#
# Reads in selected Picard metrics, generating an R-compatible TSV suitable for pre-QC analysis.
#
# There are two run modes: one pulls in timestamps using the sequencing database, one skips those two columns.
#   To load info from the sequencing database, use IntelliJ to connect to the database, thereby downloading its Oracle connection jar into $STING_HOME.  Then run:
#   /humgen/gsa-hpprojects/software/bin/jython2.5.2/jython \
#     -J-classpath $STING_HOME/dist/sam-1.48.889.jar:$STING_HOME/dist/picard-1.48.889.jar:$STING_HOME/dist/picard-private-parts-1954.jar:$STING_HOME/ojdbc6-11.2.0.1.0.jar \
#     $STING_HOME/private/python/generate_per_sample_metrics.py <bam.list> true > <output_metrics_file.tsv>
#
#  To skip the sequencing database, use the following command:
#   /humgen/gsa-hpprojects/software/bin/jython2.5.2/jython \
#     -J-classpath $STING_HOME/dist/sam-1.48.889.jar:$STING_HOME/dist/picard-1.48.889.jar:$STING_HOME/dist/picard-private-parts-1954.jar \
#     $STING_HOME/private/python/generate_per_sample_metrics.py <bam.list> > <output_metrics_file.tsv>
#
# To add a new metric:
#   - If the metric file is new to Picard, add the relevant parser to the picard-private jar
#     (see http://www.broadinstitute.org/gsa/wiki/index.php/Adding_and_updating_dependencies for details).
#   - Add the field name to the header array.
#   - Add the field data to the statement printing the data array.
#
from java.lang import *
from java.io import File,FileReader

from edu.mit.broad.picard.genotype.concordance import DbSnpMatchMetrics
from net.sf.picard.analysis import AlignmentSummaryMetrics,InsertSizeMetrics
from net.sf.picard.analysis.directed import HsMetrics
from net.sf.picard.io import IoUtil
from net.sf.picard.metrics import MetricsFile
from net.sf.picard.util import TabbedTextFileWithHeaderParser

import generate_preqc_database

import os,string,sys

def median(l):
    return sorted(l)[(len(l)+1)/2]
def mean(l):
    return float(sum(l))/len(l)

def get_all_metrics(filename):
    if not os.path.exists(filename):
        return None
    file_reader = FileReader(filename)
    metrics_file = MetricsFile()
    metrics_file.read(file_reader)
    metrics = metrics_file.getMetrics()
    file_reader.close()
    return metrics

def get_sample_summary_metrics_fields(type):
    return [field.getName() for field in type.getFields() if not field.getName().startswith('__')]

def get_sample_summary_metrics(filename,filter):
    if not os.path.exists(filename):
        return None
    file_reader = FileReader(filename)
    metrics_file = MetricsFile()
    metrics_file.read(file_reader)
    raw_metrics = metrics_file.getMetrics()
    file_reader.close()
    sampled_metrics = []
    for metric in raw_metrics:
        if filter != None:
            key,value = filter.split('=')
            if hasattr(metric,key) and getattr(metric,key).toString() == value:
                sampled_metrics.append(metric)
        else:
            sampled_metrics.append(metric)
    if len(sampled_metrics) > 1:
        raise Exception("Too many metrics to return from filename %s"%filename)
    if len(sampled_metrics) < 1:
        return None
    return sampled_metrics[0]

sample_summary_metrics_types = [ (HsMetrics,'hybrid_selection_metrics',None,None),
                                 (AlignmentSummaryMetrics,'alignment_summary_metrics','CATEGORY=PAIR',None),
                                 (InsertSizeMetrics,'insert_size_metrics','PAIR_ORIENTATION=FR','FR'),
                                 (InsertSizeMetrics,'insert_size_metrics','PAIR_ORIENTATION=RF','RF'),
                                 (InsertSizeMetrics,'insert_size_metrics','PAIR_ORIENTATION=TANDEM','TANDEM'),
                                 (DbSnpMatchMetrics,'dbsnp_matches',None,None) ]

def get_full_metrics_fields():
    sample_summary_metrics_fields = []
    # Add the initiative name from analysis_files.txt
    headers = ['INITIATIVE']
    # Add in the fingerprint lods.
    headers += ['FINGERPRINT_LODS','HAPLOTYPES_CONFIDENTLY_MATCHING']
    for metric_type in sample_summary_metrics_types:
        metric_class = metric_type[0]
        metric_suffix = metric_type[3]
        metrics_fields = get_sample_summary_metrics_fields(metric_class)
        if metric_suffix != None:
            metrics_fields = [field_name+ '_'+metric_suffix for field_name in metrics_fields]
        sample_summary_metrics_fields.extend(metrics_fields)
    return headers + sample_summary_metrics_fields

def get_full_metrics(sample,basepath):
    data = []

    # Load in the initiative from analysis_files.txt. I believe this data is lane-level, so we grab only the first row for the initiative data.
    initiative = 'NA'
    analysis_file_reader = TabbedTextFileWithHeaderParser(File(os.path.dirname(basepath)+'/analysis_files.txt'))
    for row in analysis_file_reader:
        initiative = '"%s"'%row.getField('INITIATIVE')
        break

    data += [initiative]

    fingerprinting_summary_metrics = get_all_metrics('%s.%s' % (basepath,'fingerprinting_summary_metrics'))
    
    if fingerprinting_summary_metrics != None:
        haplotypes_confidently_matching = [str(metric.HAPLOTYPES_CONFIDENTLY_MATCHING) for metric in fingerprinting_summary_metrics]
        fingerprint_lods = [str(metric.LOD_EXPECTED_SAMPLE) for metric in fingerprinting_summary_metrics]
    else:
        haplotypes_confidently_matching = []
        fingerprint_lods = []
        
    data += ['c('+string.join(fingerprint_lods,',')+')','c('+string.join(haplotypes_confidently_matching,',')+')']

    for metrics_type,metrics_extension,metrics_filter,metrics_suffix in sample_summary_metrics_types:
        metrics_pathname = '%s.%s' % (basepath,metrics_extension)
        metrics = get_sample_summary_metrics(metrics_pathname,metrics_filter)
        if metrics != None:
            data.extend([str(getattr(metrics, metrics_field_name)) for metrics_field_name in get_sample_summary_metrics_fields(metrics_type)])
        else:
            data.extend(['NA' for metrics_field_name in get_sample_summary_metrics_fields(metrics_type)])
    return data    

def main():
    include_sequence_date = False
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print >> sys.stderr, 'USAGE: %s <bam files.list> {include sequence date}'
        sys.exit(1)
    if not os.path.exists(sys.argv[1]):
        print >> sys.stderr, 'BAM list %s not found' % sys.argv[1]
        sys.exit(1)
    if len(sys.argv) == 3:
        if sys.argv[2].lower() == 'true':
            include_sequence_date = True

    bam_list_filename = sys.argv[1]

    header = ['sample']
    header.extend(get_full_metrics_fields())
    if include_sequence_date:
        header.extend(['Last_Sequenced_WR','Last_Sequenced_WR_Created_Date'])
    print string.join(header,'\t')

    # get a representative BAM file for each sample, to use as a base path.  Note that this assumes every sample corresponds to the same base path.
    bam_list = open(bam_list_filename,'r')
    samples = dict()

    for bam_filename in bam_list:
        bam_filename = bam_filename.strip()
        if bam_filename == '':
            continue
        bam_filename_tokens = bam_filename.split('/')
        project_id = bam_filename_tokens[-4]
        sample_id = bam_filename_tokens[-3]
        # Temporary hack to get through today: we currently do NOT have a way to go from a BAM list back to a sample id.
        if sample_id == 'C_III-3_A':
            sample_id = 'C III-3_A'
        samples[sample_id] = bam_filename
    bam_list.close()

    for sample_id,filename in samples.items():
        basepath = filename[:filename.rindex('.bam')]
        # Be certain to quote the metrics to ensure that spaces are properly handled by R.
        metrics = ['"%s"'%sample_id]
        metrics += get_full_metrics(sample_id,basepath)
        if include_sequence_date:
            metrics += generate_preqc_database.load_dates_from_database(project_id,sample_id)
        print >> sys.stderr, 'Processing',filename
        print string.join(metrics,'\t')

if __name__ == "__main__":
    main()

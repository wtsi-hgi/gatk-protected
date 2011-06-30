#
# Reads in selected Picard metrics, generating an R-compatible TSV suitable for pre-QC analysis.
#
# To run:
#   /humgen/gsa-hpprojects/software/bin/jython2.5.2/jython \
#     -J-classpath $STING_HOME/dist/sam-1.47.869.jar:$STING_HOME/dist/picard-1.47.869.jar:$STING_HOME/dist/picard-private-parts-1941.jar \
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
from net.sf.picard.metrics import MetricsFile

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
        raise Exception("Too few metrics to return from filename %s"%filename)
    return sampled_metrics[0]

sample_summary_metrics_types = [ (HsMetrics,'hybrid_selection_metrics',None),
                                 (AlignmentSummaryMetrics,'alignment_summary_metrics','CATEGORY=PAIR'),
                                 (InsertSizeMetrics, 'insert_size_metrics',None),
                                 (DbSnpMatchMetrics, 'dbsnp_matches',None) ]

def get_full_metrics_fields():
    headers = ['FINGERPRINT_LODS','HAPLOTYPES_CONFIDENTLY_MATCHING']
    for metric_type in sample_summary_metrics_types:
        headers.extend(get_sample_summary_metrics_fields(metric_type[0]))
    return headers

def get_full_metrics(sample,basepath):
    fingerprinting_summary_metrics = get_all_metrics('%s.%s' % (basepath,'fingerprinting_summary_metrics'))
    
    if fingerprinting_summary_metrics != None:
        haplotypes_confidently_matching = [str(metric.HAPLOTYPES_CONFIDENTLY_MATCHING) for metric in fingerprinting_summary_metrics]
        fingerprint_lods = [str(metric.LOD_EXPECTED_SAMPLE) for metric in fingerprinting_summary_metrics]
    else:
        haplotypes_confidently_matching = []
        fingerprint_lods = []
        
    data = ['c('+string.join(fingerprint_lods,',')+')','c('+string.join(haplotypes_confidently_matching,',')+')']

    for metrics_type,metrics_extension,metrics_filter in sample_summary_metrics_types:
        metrics_pathname = '%s.%s' % (basepath,metrics_extension)
        if os.path.exists(metrics_pathname):
            metrics = get_sample_summary_metrics(metrics_pathname,metrics_filter)
            data.extend([str(getattr(metrics, metrics_field_name)) for metrics_field_name in get_sample_summary_metrics_fields(metrics_type)])
        else:
            data.extend(['NA' for metrics_field_name in get_sample_summary_metrics_fields(metrics_type)])
    return data    

def main():
    if len(sys.argv) != 2:
        print 'USAGE: %s <bam files.list>'
        sys.exit(1)
    if not os.path.exists(sys.argv[1]):
        print 'BAM list %s not found' % sys.argv[1]
        sys.exit(1)

    bam_list_filename = sys.argv[1]

    header = ['sample']
    for metric_type in sample_summary_metrics_types:
        header.extend(get_sample_summary_metrics_fields(metric_type[0]))
    print string.join(header,'\t')

    # get a representative BAM file for each sample, to use as a base path.  Note that this assumes every sample corresponds to the same base path.
    bam_list = open(bam_list_filename,'r')
    samples = dict()

    for bam_filename in bam_list:
        bam_filename = bam_filename.strip()
        if bam_filename == '':
            continue
        bam_filename_tokens = bam_filename.split('/')
        sample_id = bam_filename_tokens[len(bam_filename_tokens)-3]
        samples[sample_id] = bam_filename
    bam_list.close()

    for sample_id,filename in samples.items():
        basepath = filename[:filename.rindex('.bam')]
        print string.join([sample_id]+get_full_metrics(sample_id,basepath),'\t')

if __name__ == "__main__":
    main()

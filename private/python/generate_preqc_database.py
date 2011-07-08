#
# To run:
#   /humgen/gsa-hpprojects/software/bin/jython2.5.2/jython \
#     -J-classpath $STING_HOME/lib/poi-3.8-beta3.jar:$STING_HOME/lib/poi-ooxml-3.8-beta3.jar:$STING_HOME/lib/poi-ooxml-schemas-3.8-beta3.jar:$STING_HOME/lib/xmlbeans-2.3.0.jar:$STING_HOME/lib/dom4j-1.6.1.jar:$STING_HOME/lib/sam-1.47.869.jar:$STING_HOME/lib/picard-1.47.869.jar:$STING_HOME/lib/picard-private-parts-1941.jar:$STING_HOME/dist/GenomeAnalysisTK.jar:$STING_HOME/ojdbc6-11.2.0.1.0.jar \
#     generate_preqc_database.py <input file>
#
import os,string,sys

from parse_pm_input import project_file_reader
from generate_per_sample_metrics import get_full_metrics_fields,get_full_metrics

from java.io import File
from java.lang import Class,Exception
from java.sql import *
from java.text import SimpleDateFormat

functional_classes = ['all','missense','nonsense','silent']
novelties = ['all','known','novel']

count_variants_columns = ['CompRod','EvalRod','FunctionalClass','Novelty','nProcessedLoci','nCalledLoci','nRefLoci','nVariantLoci','variantRate','variantRatePerBp',
                          'nSNPs','nMNPs','nInsertions','nDeletions','nComplex','nNoCalls','nHets','nHomRef','nHomVar','nSingletons','nHomDerived','heterozygosity',
                          'heterozygosityPerBp','hetHomRatio','indelRate','indelRatePerBp','deletionInsertionRatio']
titv_variant_evaluator_columns = ['nTi','nTv','tiTvRatio','nTiInComp','nTvInComp','TiTvRatioStandard','nTiDerived','nTvDerived','tiTvDerivedRatio']

def load_dates_from_database(project,sample):
    Class.forName("oracle.jdbc.OracleDriver")
    url = "jdbc:oracle:thin:REPORTING/REPORTING@//ora01:1521/SEQPROD"
    con = DriverManager.getConnection(url)
    stmt = con.createStatement()
    sql = 'select "Last Sequenced WR","Last Sequenced WR Created Date" from ILLUMINA_SAMPLE_STATUS_AGG where "Project" = \'%s\' and "Sample" = \'%s\'' % (project,sample)
    rs = stmt.executeQuery(sql)

    last_sequenced_wr = None
    last_sequenced_wr_created_date = None

    has_results = rs.next()

    if has_results and rs.getString(1) != None:
        last_sequenced_wr = rs.getString(1)
    if has_results and rs.getString(2) != None:
        date_formatter = SimpleDateFormat('yyyy-MM-dd')
        last_sequenced_wr_created_date = date_formatter.format(rs.getDate(2))

    rs.close()
    stmt.close()
    con.close()
    return last_sequenced_wr,last_sequenced_wr_created_date

def generate_project_files_from_filtered_annotated_vcfs(vcf_list_file):
    vcf_list = open(vcf_list_file,'r')
    for filename in vcf_list:
        filename = filename.strip()
        project_dirname = os.path.dirname(filename)
        project_basename = os.path.basename(filename)[:os.path.basename(filename).find('.')]
        project_filename = '%s/%s.tsv'%(project_dirname,project_basename)
        if not os.path.exists(project_filename):
            print >> sys.stderr,'WARNING: Unable to find file',project_filename
            continue
        try:
            for squid,sample,latest_version in project_file_reader(project_filename):
                yield project_basename,squid,sample,latest_version
        except: 
            print >> sys.stderr,'WARNING: Cannot decode file format for',project_filename
    vcf_list.close()

def main():
    from org.broadinstitute.sting.gatk.report import GATKReportParser

    if len(sys.argv) != 3:
        print >> sys.stderr, 'USAGE: %s <vcf.list> <variant eval outputs dir>'
        sys.exit(1)
    if not os.path.exists(sys.argv[1]):
        print >> sys.stderr, 'VCF list file %s not found' % sys.argv[1]
        sys.exit(1)
    if not os.path.exists(sys.argv[2]):
        print >> sys.stderr, 'VariantEval output directory %s not found' % sys.argv[2]
        sys.exit(1)

    vcf_list = sys.argv[1]
    variant_eval_base = sys.argv[2]

    # print out headers
    print string.join(['project','squid','sample']+count_variants_columns+titv_variant_evaluator_columns+get_full_metrics_fields()+['Last_Sequenced_WR','Last_Sequenced_WR_Created_Date'],'\t')

    for project,squid,sample,latest_version in generate_project_files_from_filtered_annotated_vcfs(vcf_list):
        print >> sys.stderr, 'processing project = %s, squid id = %s, sample = %s'%(project,squid,sample)
        variant_eval_filename = '%s/%s.cleaned.snps_and_indels.filtered.annotated.perSample.exons.eval'%(variant_eval_base,project)
        if not os.path.exists(variant_eval_filename):
            print 'WARNING: varianteval path %s does not exist'%variant_eval_filename
            continue
        report_parser = GATKReportParser()
        report_parser.parse(File(variant_eval_filename))

        last_sequenced_wr,last_sequenced_wr_created_date = load_dates_from_database(squid,sample)

        for functional_class in functional_classes:
            for novelty in novelties:
                columns = [project,squid,sample]
                columns.extend([report_parser.getValue('CountVariants',['dbsnp','eval',functional_class,novelty,sample],column) for column in count_variants_columns])
                columns.extend([report_parser.getValue('TiTvVariantEvaluator',['dbsnp','eval',functional_class,novelty,sample],column) for column in titv_variant_evaluator_columns])
                columns.extend(get_full_metrics(sample,'/seq/picard_aggregation/%s/%s/v%s/%s'%(squid,sample,latest_version,sample)))
                columns.extend([last_sequenced_wr,last_sequenced_wr_created_date])
                # replace any Nones in the columns with the text 'NA'
                columns = [column if column != None else 'NA' for column in columns]
                print string.join(columns,'\t')
            
if __name__ == "__main__":
    main()

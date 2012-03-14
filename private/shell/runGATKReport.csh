#!/bin/tcsh

source /broad/tools/scripts/useuse

reuse Python-2.5

setenv DIR /humgen/gsa-hpprojects/GATK/reports 
setenv ARCHIVE_DIR $DIR/archive
setenv DATE `date +"%m_%d_%Y"`
setenv ARCHIVE $ARCHIVE_DIR/$DATE
setenv GATK ~/dev/GenomeAnalysisTK/unstable/private
setenv REPORT_TXT $DIR/report.txt

rm -f $REPORT_TXT 

cd $DIR

echo "\n####################\nArchiving recently submitted jobs" >> $REPORT_TXT
python $GATK/python/analyzeRunReports.py archive $DIR/submitted -o $ARCHIVE.gz -D >> $REPORT_TXT

echo "\n####################\nLoading to DB" >> $REPORT_TXT
python $GATK/python/analyzeRunReports.py loadToDB $ARCHIVE.gz >> $REPORT_TXT

cat $REPORT_TXT



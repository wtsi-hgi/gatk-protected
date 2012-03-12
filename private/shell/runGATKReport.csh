#!/bin/tcsh

source /broad/tools/scripts/useuse

reuse Python-2.5
use R-2.11

setenv DIR /humgen/gsa-hpprojects/GATK/reports 
setenv ARCHIVE_DIR $DIR/archive
setenv SUMMARY_DIR $DIR/summaries
setenv DATE `date +"%m_%d_%Y"`
setenv ARCHIVE $ARCHIVE_DIR/$DATE
setenv SUMMARY $SUMMARY_DIR/$DATE
setenv GATK ~/dev/GenomeAnalysisTK/unstable/private
setenv GATK_RELEASE_VERSION `ls -l /humgen/gsa-hpprojects/GATK/bin/current | sed 's/.*GenomeAnalysisTK-\([0-9]*\.[0-9]*-\).*/\1/'`
setenv REPORT_TXT $DIR/report.txt

rm -f $REPORT_TXT 

cd $DIR

echo "\n####################\nArchiving recently submitted jobs" >> $REPORT_TXT
python $GATK/python/analyzeRunReports.py archive $DIR/submitted -o $ARCHIVE.gz -D >> $REPORT_TXT

# echo "\n####################\nReleased version ($GATK_RELEASE_VERSION), all runs" >> $REPORT_TXT
# python $GATK/python/analyzeRunReports.py summary $ARCHIVE_DIR/*.gz --rev $GATK_RELEASE_VERSION >> $REPORT_TXT
# python $GATK/python/analyzeRunReports.py exceptions $ARCHIVE_DIR/*.gz -E sting --rev $GATK_RELEASE_VERSION >> $REPORT_TXT

echo "\n####################\nLoading to DB" >> $REPORT_TXT
python $GATK/python/analyzeRunReports.py loadToDB $ARCHIVE.gz >> $REPORT_TXT

echo "\n####################\nLast day, all versions summary" >> $REPORT_TXT
python $GATK/python/analyzeRunReports.py summary $ARCHIVE.gz --no-dev >> $REPORT_TXT

echo "\n####################\nLast day exceptions for rev $GATK_RELEASE_VERSION" >> $REPORT_TXT
python $GATK/python/analyzeRunReports.py exceptions $ARCHIVE.gz -E sting --no-dev --rev $GATK_RELEASE_VERSION >> $REPORT_TXT

#echo "GATK daily run report" | mutt -a $SUMMARY.30_days.pdf -a $SUMMARY.360_days.pdf -a $SUMMARY.7_days.pdf -s "GATK Run report PDFs for $DATE" gsamembers
#cat $REPORT_TXT | mutt -a $REPORT_TXT -a $SUMMARY.30_days.pdf -a $SUMMARY.360_days.pdf -s "GATK run report for $DATE" gsamembers
#cat $REPORT_TXT | mutt -a $REPORT_TXT -s "GATK run report for $DATE" gsamembers
cat $REPORT_TXT



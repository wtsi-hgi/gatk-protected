#!/bin/tcsh

reuse Python-2.6

setenv DATE `date +"%m_%d_%Y"`
setenv ROOT /humgen/gsa-hpprojects/GATK/reports/s3
setenv GATK /local/gsa-engineering/cron_clones/unstable
setenv DOWNLOAD_ROOT /local/dev/GATKLogs
setenv DOWNLOAD_DIR $DOWNLOAD_ROOT/$DATE
setenv DIR $ROOT/archive/$DATE

#
# manageGATKS3Logs.py copied from GATK repository, do not modify it here
#
setenv BASE "python $GATK/private/GATKLogs/manageGATKS3Logs.py -b GATK_Run_Reports -s $ROOT/s3cmd-1.0.0/s3cmd -d $DOWNLOAD_DIR -p 10 -g 100"
setenv LSFILE $DOWNLOAD_ROOT/files_$DATE.ls

$BASE ls $LSFILE
$BASE move $LSFILE progress_$DATE.log
echo 'Done:', `date`

echo "\n####################\nCreating raw tgz"
tar -czf $DIR.raw_files.gz -C $DOWNLOAD_ROOT $DATE
echo 'Done:', `date`

echo "\n####################\nArchiving"
python $GATK/private/GATKLogs/analyzeRunReports.py archive $DOWNLOAD_DIR -o $DIR.gz -D
echo 'Done:', `date`

echo "\n####################\nLoading to DB"
python $GATK/private/GATKLogs/analyzeRunReports.py loadToDB $DIR.gz
echo 'Done:', `date`

# if the dir is empty we proceed
rmdir --ignore-fail-on-non-empty $DIR
rm -f progress_$DATE.log


#!/bin/tcsh

foreach mode (setupDB loadToDB)
python private/python/analyzeRunReports.py $mode --updateFreq 10000 /humgen/gsa-hpprojects/GATK/reports/archive/*.gz /humgen/gsa-hpprojects/GATK/reports/s3/archive/*.gz
end

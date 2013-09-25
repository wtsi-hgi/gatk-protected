#!/bin/tcsh

use Java-1.7

set args = ""

if ( $1 == "local" ) then
set args = " -useLocal"
shift
endif

if ( $1 == "dev" ) then
set args = "$args -dbToUse DEV"
shift
endif

set logDirectory = /humgen/gsa-hpprojects/NA12878Collection/knowledgeBase/logs
set reviewDirectory = /humgen/gsa-hpprojects/NA12878Collection/knowledgeBase/reviewsBackup
set DATE = `date +"%m_%d_%Y"`
set reviewVCFName = $reviewDirectory/$DATE.reviews.vcf
set BUNDLE = /humgen/gsa-hpprojects/GATK/bundle/current

java -Xmx2g -jar ${1} -R $BUNDLE/b37/human_g1k_v37.fasta $args -T NA12878KnowledgeBaseServer -reviewsFile $reviewVCFName -maxRuntime 1430 -maxRuntimeUnits MINUTES -dbToUse PRODUCTION -log $logDirectory/$DATE.log

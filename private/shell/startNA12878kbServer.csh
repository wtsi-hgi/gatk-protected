#!/bin/tcsh

set args = ""

if ( $1 == "local" ) then
set args = " -useLocal"
shift
endif

if ( $1 == "dev" ) then
set args = "$args -dbToUse DEV"
shift
endif

set reviewDirectory = /humgen/gsa-hpprojects/NA12878Collection/knowledgeBase/reviewsBackup
set DATE = `date +"%m_%d_%Y"`
set reviewVCFName = $reviewDirectory/$DATE.reviews.vcf

java -Xmx2g -jar dist/GenomeAnalysisTK.jar -R $BUNDLE/b37/human_g1k_v37.fasta $args -T NA12878KnowledgeBaseServer -reviewsFile $reviewVCFName -maxRuntime 1430 -maxRuntimeUnits MINUTES -dbToUse PRODUCTION

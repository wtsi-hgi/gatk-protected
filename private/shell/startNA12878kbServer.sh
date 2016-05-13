#!/bin/bash

source /broad/software/scripts/useuse
use Java-1.8

args=""

if [ "$1" == "local" ] 
then
    args=" -useLocal"
    shift
fi

if [ "$1" == "dev" ] 
then
    args="${args} -dbToUse DEV"
    shift
fi

logDirectory="/humgen/gsa-hpprojects/NA12878Collection/knowledgeBase/logs"
reviewDirectory="/humgen/gsa-hpprojects/NA12878Collection/knowledgeBase/reviewsBackup"
DATE=`date +"%m_%d_%Y"`
reviewVCFName="${reviewDirectory}/${DATE}.reviews.vcf"
BUNDLE="/humgen/gsa-hpprojects/GATK/bundle/current"

java -Xmx2g -jar ${1} -R ${BUNDLE}/b37/human_g1k_v37.fasta $args -T NA12878KnowledgeBaseServer -reviewsFile ${reviewVCFName} -maxRuntime 1430 -maxRuntimeUnits MINUTES -dbToUse PRODUCTION -log ${logDirectory}/${DATE}.log

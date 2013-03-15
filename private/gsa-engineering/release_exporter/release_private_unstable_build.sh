#!/bin/bash
#
# Publish a full GATK jar that includes private to an internal directory
# for GSA group use.
#

if [ $# -ne 2 ]
then
    echo "Usage: $0 ivy_cache_directory temp_directory"
    exit 1
fi

IVY_CACHE_DIR="$1"
TEMP_DIR="$2"
DESTINATION_DIR="/humgen/gsa-hpprojects/GATK/private_unstable_builds"
GATK_VERSION=`git describe --long`
PACKAGE_OUTPUT_DIR="dist/packages/GenomeAnalysisTK-${GATK_VERSION}"
PACKAGE_JAR_RELATIVE_PATH="GenomeAnalysisTK-${GATK_VERSION}/GenomeAnalysisTK.jar"
CURRENT_UNSTABLE_JAR_LINK_NAME="GenomeAnalysisTK_latest_unstable.jar"

ant clean package.gatk.all "-Djava.io.tmpdir=${TEMP_DIR}" "-Divy.home=${IVY_CACHE_DIR}"

if [ $? -ne 0 ]
then
    echo "$0: failed to package GATK jar"
    exit 1
fi

cp -r "${PACKAGE_OUTPUT_DIR}" "${DESTINATION_DIR}" && \
cd "${DESTINATION_DIR}" && \
ln -fs "${PACKAGE_JAR_RELATIVE_PATH}" "${CURRENT_UNSTABLE_JAR_LINK_NAME}"

if [ $? -ne 0 ]
then
    echo "$0: failed to publish private GATK jar for version ${GATK_VERSION} to directory ${DESTINATION_DIR}"
    exit 1
fi

echo "$0: successfully published private GATK jar for version ${GATK_VERSION} to directory ${DESTINATION_DIR}"
exit 0

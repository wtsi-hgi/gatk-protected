#!/bin/bash
#
# Publish full GATK and Queue jars (including private classes) to an internal directory
# for GSA group use.
#

WORKING_DIR=`pwd`
TEMP_DIR="tmp"
BUILD_VERSION=`git describe --long`
GATK_DESTINATION_ROOT="/humgen/gsa-hpprojects/GATK/private_unstable_builds"
QUEUE_DESTINATION_ROOT="/humgen/gsa-hpprojects/Queue/private_unstable_builds"
GATK_DESTINATION_SUBDIR="GenomeAnalysisTK-${BUILD_VERSION}"
QUEUE_DESTINATION_SUBDIR="Queue-${BUILD_VERSION}"
GATK_PACKAGE_JAR="target/GenomeAnalysisTK.jar"
QUEUE_PACKAGE_JAR="target/Queue.jar"
CURRENT_UNSTABLE_GATK_JAR_LINK_NAME="GenomeAnalysisTK_latest_unstable.jar"
CURRENT_UNSTABLE_QUEUE_JAR_LINK_NAME="Queue_latest_unstable.jar"

mkdir "${TEMP_DIR}"

mvn clean && mvn package "-Djava.io.tmpdir=${TEMP_DIR}"

if [ $? -ne 0 ]
then
    echo "$0: failed to package GATK/Queue jars"
    exit 1
fi

mkdir "${GATK_DESTINATION_ROOT}/${GATK_DESTINATION_SUBDIR}" && \
cp "${GATK_PACKAGE_JAR}" "${GATK_DESTINATION_ROOT}/${GATK_DESTINATION_SUBDIR}/GenomeAnalysisTK.jar" && \
cd "${GATK_DESTINATION_ROOT}" && \
ln -fs "${GATK_DESTINATION_SUBDIR}/GenomeAnalysisTK.jar" "${CURRENT_UNSTABLE_GATK_JAR_LINK_NAME}"

if [ $? -ne 0 ]
then
    echo "$0: failed to publish private GATK jar for version ${BUILD_VERSION} to directory ${GATK_DESTINATION_ROOT}"
    exit 1
fi

cd "${WORKING_DIR}"

mkdir "${QUEUE_DESTINATION_ROOT}/${QUEUE_DESTINATION_SUBDIR}" && \
cp "${QUEUE_PACKAGE_JAR}" "${QUEUE_DESTINATION_ROOT}/${QUEUE_DESTINATION_SUBDIR}/Queue.jar" && \
cd "${QUEUE_DESTINATION_ROOT}" && \
ln -fs "${QUEUE_DESTINATION_SUBDIR}/Queue.jar" "${CURRENT_UNSTABLE_QUEUE_JAR_LINK_NAME}"

if [ $? -ne 0 ]
then
    echo "$0: failed to publish private Queue jar for version ${BUILD_VERSION} to directory ${QUEUE_DESTINATION_ROOT}"
    exit 1
fi

echo "$0: successfully published private GATK and Queue jars for version ${BUILD_VERSION} to directories ${GATK_DESTINATION_ROOT} and ${QUEUE_DESTINATION_ROOT}"
exit 0

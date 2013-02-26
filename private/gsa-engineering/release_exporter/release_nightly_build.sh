#!/bin/bash

if [ $# -ne 2 ]
then
    echo "Usage: $0 ivy_cache_directory temp_directory"
    exit 1
fi

IVY_CACHE_DIR="$1"
TEMP_DIR="$2"
DESTINATION_DIR="/humgen/gsa-hpprojects/GATK/nightly_builds"
PACKAGE_OUTPUT_DIR="dist/packages"
GATKDOCS_OUTPUT_DIR="gatkdocs"

TIMESTAMP=`date '+%Y-%m-%d'`
HASH_PREFIX=`git describe --long | awk -F'-' '{ print $3; }'`
NIGHTLY_BUILD_VERSION="nightly-${TIMESTAMP}-${HASH_PREFIX}"

GATK_ARCHIVE_NAME="GenomeAnalysisTK-${NIGHTLY_BUILD_VERSION}.tar.bz2"
GATKDOCS_ARCHIVE_NAME="Gatkdocs-${NIGHTLY_BUILD_VERSION}.tar.bz2"
CURRENT_NIGHTLY_GATK_LINK_NAME="nightly_gatk.tar.bz2"
CURRENT_NIGHTLY_GATKDOCS_LINK_NAME="nightly_gatkdocs.tar.bz2"

ant clean package.gatk.full "-Djava.io.tmpdir=${TEMP_DIR}" \
                            "-Divy.home=${IVY_CACHE_DIR}" \
                            "-Dbuild.version=${NIGHTLY_BUILD_VERSION}"

if [ $? -ne 0 ]
then
    echo "Failed to package GATK jar for nightly build ${NIGHTLY_BUILD_VERSION}"
    exit 1
fi

ant gatkdocs "-Djava.io.tmpdir=${TEMP_DIR}" \
             "-Divy.home=${IVY_CACHE_DIR}" \
             "-Dbuild.version=${NIGHTLY_BUILD_VERSION}"

if [ $? -ne 0 ]
then
    echo "Failed generate gatkdocs for nightly build ${NIGHTLY_BUILD_VERSION}"
    exit 1
fi

tar -c -j -f "${GATKDOCS_ARCHIVE_NAME}" "${GATKDOCS_OUTPUT_DIR}"

if [ $? -ne 0 ]
then
    echo "Failed create gatkdocs compressed archive for nightly build ${NIGHTLY_BUILD_VERSION}"
    exit 1
fi

cp "${PACKAGE_OUTPUT_DIR}/${GATK_ARCHIVE_NAME}" "${DESTINATION_DIR}" && \
cp "${GATKDOCS_ARCHIVE_NAME}" "${DESTINATION_DIR}" && \
cd "${DESTINATION_DIR}" && \
ln -fs "${GATK_ARCHIVE_NAME}" "${CURRENT_NIGHTLY_GATK_LINK_NAME}" && \
ln -fs "${GATKDOCS_ARCHIVE_NAME}" "${CURRENT_NIGHTLY_GATKDOCS_LINK_NAME}"

if [ $? -ne 0 ]
then
    echo "Failed to publish nightly build ${NIGHTLY_BUILD_VERSION} to directory ${DESTINATION_DIR}"
    exit 1
fi

echo "Successfully published GATK nightly build version ${NIGHTLY_BUILD_VERSION}"
exit 0
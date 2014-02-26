#!/bin/sh

ORIGINAL_PATH=$1
ORIGINAL_NAME=`basename $ORIGINAL_PATH`
ORIGINAL_PATTERN='([A-Za-z]+)-(.+).tar.bz2'

if [[ "$ORIGINAL_NAME" =~ $ORIGINAL_PATTERN ]]; then
    PACKAGE_NAME=${BASH_REMATCH[1]}
    GIT_VERSION=${BASH_REMATCH[2]}
    ORIGINAL_ROOT=${PACKAGE_NAME}-${GIT_VERSION}
else
    echo "Unable to regex package name and git version: $1" >&2
    exit 1
fi

if [ "$PACKAGE_NAME" == "GenomeAnalysisTK" ]; then
    PROJECT_NAME=GATK
elif [ "$PACKAGE_NAME" == "Queue" ]; then
    PROJECT_NAME=Queue
else
    echo "Unknown project $PACKAGE_NAME" >&2
fi

PROJECT_DIR=/humgen/gsa-hpprojects/${PROJECT_NAME}/bin
EXPANDED_DIR=${PROJECT_DIR}/${ORIGINAL_ROOT}
CURRENT_DIR_LINK=${PROJECT_DIR}/current
LATEST_LINK=${PROJECT_DIR}/${PACKAGE_NAME}-latest.tar.bz2

set -e -x

cp ${ORIGINAL_PATH} ${PROJECT_DIR}/.
mkdir ${EXPANDED_DIR}
tar -xf ${ORIGINAL_PATH} -C ${EXPANDED_DIR}
ln -sf ${EXPANDED_DIR} ${CURRENT_DIR_LINK}
ln -sf ${PROJECT_DIR}/${ORIGINAL_NAME} ${LATEST_LINK}

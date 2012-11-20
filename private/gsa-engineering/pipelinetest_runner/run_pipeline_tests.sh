#!/bin/bash

TEST_NAME="$1"
IVY_CACHE_DIR="$2"
TEMP_DIR="$3"
PIPELINE_TEST_ROOT_DIR=/humgen/gsa-scr1/gsa-engineering/pipeline_test_working_directory
PIPELINE_TEST_WORKING_DIR="${PIPELINE_TEST_ROOT_DIR}/${TEST_NAME}"
BAMBOO_CLONE=`pwd`

if [ $# -ne 3 ]
then
    echo "$0: Usage: $0 test_name ivy_cache_dir temp_dir"
    exit 1
fi

if [ "${TEST_NAME}" != "stable" -a "${TEST_NAME}" != "unstable" ]
then
    echo "$0: Invalid test name: ${TEST_NAME}"
    exit 1
fi

if [ ! -d "${PIPELINE_TEST_ROOT_DIR}" ]
then
    echo "$0: ${PIPELINE_TEST_ROOT_DIR} does not exist. Please create it as user gsa-engineering."
    exit 1
fi

if [ ! -d "${IVY_CACHE_DIR}" ]
then
    echo "$0: Ivy cache directory ${IVY_CACHE_DIR} does not exist. Bamboo should have created this directory."
    exit 1
fi

if [ ! -d "${TEMP_DIR}" ]
then
    echo "$0: Temp directory ${TEMP_DIR} does not exist. Please create it."
    exit 1
fi

if [ -d "${PIPELINE_TEST_WORKING_DIR}" ]
then
    rm -rf "${PIPELINE_TEST_WORKING_DIR}"
fi

echo "$0: Cloning ${BAMBOO_CLONE} into ${PIPELINE_TEST_WORKING_DIR}"
git clone "${BAMBOO_CLONE}" "${PIPELINE_TEST_WORKING_DIR}"

if [ $? -ne 0 ]
then
    echo "$0: Failed to clone ${BAMBOO_CLONE} into ${PIPELINE_TEST_WORKING_DIR}"
    exit 1
fi

cd "${PIPELINE_TEST_WORKING_DIR}"
ant clean pipelinetestrun -Djava.io.tmpdir=${TEMP_DIR} -Divy.home=${IVY_CACHE_DIR} -Dhalt=no

if [ $? -ne 0 ]
then
    echo "$0: ${TEST_NAME} pipelinetestrun failed."
    exit 1
fi

echo "$0: ${TEST_NAME} pipelinetestrun succeeded."

mkdir "${BAMBOO_CLONE}/build" && cp -r "${PIPELINE_TEST_WORKING_DIR}/build/report" "${BAMBOO_CLONE}/build"

if [ $? -ne 0 ]
then
    echo "$0: Failed to copy test results from ${PIPELINE_TEST_WORKING_DIR}/build/report into ${BAMBOO_CLONE}/build"
    exit 1
fi

exit 0


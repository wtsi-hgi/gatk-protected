#!/bin/bash

TEST_NAME="$1"
TEMP_DIR="$2"
QUEUE_TEST_ROOT_DIR=/humgen/gsa-scr1/gsa-engineering/queue_test_working_directory
QUEUE_TEST_WORKING_DIR="${QUEUE_TEST_ROOT_DIR}/${TEST_NAME}"
BAMBOO_CLONE=`pwd`

if [ $# -ne 2 ]
then
    echo "$0: Usage: $0 test_name temp_dir"
    exit 1
fi

if [ "${TEST_NAME}" != "stable" -a "${TEST_NAME}" != "unstable" ]
then
    echo "$0: Invalid test name: ${TEST_NAME}"
    exit 1
fi

if [ ! -d "${QUEUE_TEST_ROOT_DIR}" ]
then
    echo "$0: ${QUEUE_TEST_ROOT_DIR} does not exist. Please create it as user gsa-engineering."
    exit 1
fi

if [ ! -d "${TEMP_DIR}" ]
then
    mkdir "${TEMP_DIR}"
fi

if [ -d "${QUEUE_TEST_WORKING_DIR}" ]
then
    rm -rf "${QUEUE_TEST_WORKING_DIR}"
fi

echo "$0: Cloning ${BAMBOO_CLONE} into ${QUEUE_TEST_WORKING_DIR}"
git clone "${BAMBOO_CLONE}" "${QUEUE_TEST_WORKING_DIR}"

if [ $? -ne 0 ]
then
    echo "$0: Failed to clone ${BAMBOO_CLONE} into ${QUEUE_TEST_WORKING_DIR}"
    exit 1
fi

cd "${QUEUE_TEST_WORKING_DIR}"
git checkout -f HEAD
mvn verify -Dsting.queuetests.skipped=false -Dsting.queuetests.run=true -Djava.io.tmpdir=${TEMP_DIR}

if [ $? -ne 0 ]
then
    echo "$0: ${TEST_NAME} queue test run failed."
else
    echo "$0: ${TEST_NAME} queue test run succeeded."
fi

mkdir "${BAMBOO_CLONE}/test_results"
RESULT_COUNT=0
for test_result in `find . -name 'TEST-*.xml'`
do
    RESULT_COUNT=`expr ${RESULT_COUNT} + 1`
    RESULT_BASENAME=`basename "${test_result}" ".xml"`
    cp "${test_result}" "${BAMBOO_CLONE}/test_results/${RESULT_BASENAME}_${RESULT_COUNT}.xml"

    if [ $? -ne 0 ]
    then
        echo "$0: Failed to copy test result ${test_result} into ${BAMBOO_CLONE}/test_results"
        exit 1
    fi
done

exit 0

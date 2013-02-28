#!/bin/bash

if [ $# -ne 3 ]
then
    echo "Usage: $0 test_class_suffix job_queue temp_dir"
    echo "Example test class suffixes: UnitTest, IntegrationTest, PipelineTest"
    exit 1
fi

TEST_CLASS_SUFFIX="$1"
JOB_QUEUE="$2"
TEMP_DIR="$3"
BAMBOO_CLONE=`pwd`
# name of the bamboo build directory can serve as a unique id for this run
BAMBOO_BUILD_ID=`basename "${BAMBOO_CLONE}"`
TEST_ROOT_WORKING_DIR="/humgen/gsa-scr1/gsa-engineering/parallel_tests_working_directory/${BAMBOO_BUILD_ID}"
TEST_CLONE="${TEST_ROOT_WORKING_DIR}/test_clone"
IVY_CACHE="${TEST_ROOT_WORKING_DIR}/ivy_cache"
JOB_OUTPUT_DIR="${TEST_ROOT_WORKING_DIR}/job_output"
JOB_MEMORY="4"
# Timeout in seconds before we give up and conclude that our jobs are stuck
GLOBAL_TIMEOUT=3600
# We check job status every JOB_POLL_INTERVAL seconds
JOB_POLL_INTERVAL=30

if [ -d "${TEST_ROOT_WORKING_DIR}" ]
then
    rm -rf "${TEST_ROOT_WORKING_DIR}"
fi

mkdir "${TEST_ROOT_WORKING_DIR}"
mkdir "${IVY_CACHE}"
mkdir "${JOB_OUTPUT_DIR}"

echo "$0: Cloning ${BAMBOO_CLONE} into ${TEST_CLONE}"
git clone "${BAMBOO_CLONE}" "${TEST_CLONE}"

if [ $? -ne 0 ]
then
    echo "$0: failed to clone ${BAMBOO_CLONE} into ${TEST_CLONE}"
    exit 1
fi

cd "${TEST_CLONE}"
ant clean test.compile -Divy.home="${IVY_CACHE}" -Djava.io.tmpdir="${TEMP_DIR}"

if [ $? -ne 0 ]
then
    echo "$0: test.compile failed"
    exit 1
fi

echo "Dispatching jobs for ${BAMBOO_BUILD_ID}"
NUM_JOBS=0

for test_class_file in `find . -name "*${TEST_CLASS_SUFFIX}.java"`
do
    test_class=`basename "${test_class_file}" ".java"`

    bsub -P "${BAMBOO_BUILD_ID}" \
         -q "${JOB_QUEUE}" \
         -R "rusage[mem=${JOB_MEMORY}] select[tmp>100]" \
         -o "${JOB_OUTPUT_DIR}/${test_class}.job.out" \
         ant runtestonly -Dsingle="${test_class}" -Divy.home="${IVY_CACHE}" -Djava.io.tmpdir="${TEMP_DIR}"
    
    if [ $? -ne 0 ]
    then
        echo "Failed to dispatch job for test class ${test_class}"
        # TODO: retry in this case
        exit 1
    fi

    NUM_JOBS=`expr "${NUM_JOBS}" + 1`
done

echo "Dispatched ${NUM_JOBS} jobs to job queue ${JOB_QUEUE}"

ELAPSED_SECONDS=0
while [ "${ELAPSED_SECONDS}" -lt "${GLOBAL_TIMEOUT}" ]
do
    sleep "${JOB_POLL_INTERVAL}" 
    ELAPSED_SECONDS=`expr "${ELAPSED_SECONDS}" + "${JOB_POLL_INTERVAL}"`

    NUM_OUTSTANDING_JOBS=`bjobs -p -r -P "${BAMBOO_BUILD_ID}" 2> /dev/null | wc -l`
    if [ "${NUM_OUTSTANDING_JOBS}" -eq 0 ]
    then
        echo "$0: All jobs for ${BAMBOO_BUILD_ID} complete, saving results into ${BAMBOO_CLONE}"
        
        cp -r "build/report" "${BAMBOO_CLONE}/parallel_test_results"
        if [ $? -ne 0 ]
        then
            echo "$0: failed to copy parallel test results into ${BAMBOO_CLONE}/parallel_test_results"
            exit 1
        fi

        echo "$0: parallel test results saved"
        # TODO: delete job working directory before exiting (keep it for now for debugging purposes)
        exit 0
    fi
done

echo "$0: Timeout of ${GLOBAL_TIMEOUT} seconds reached before jobs completed, giving up"
# TODO: more graceful shutdown including bkill of outstanding jobs
exit 1


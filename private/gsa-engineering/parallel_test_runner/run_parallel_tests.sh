#!/bin/bash
#
# Simple implementation of class-level test suite parallelism. Dispatches one job
# to the farm per test class, monitors jobs for completion, and saves the aggregated
# results for display in Bamboo.
#
# Usage: run_parallel_tests.sh test_class_suffix job_queue timeout bamboo_build_number temp_dir
#
# Arguments:
#
# test_class_suffix: run tests from test classes ending in this suffix
#                    (eg., UnitTest, IntegrationTest, etc.)
#
# job_queue: LSF queue to which to submit jobs
#
# timeout: if all jobs don't complete within this many seconds, terminate all jobs
#          and exit with an error
#
# bamboo_build_number: build number of the Bamboo build that is invoking this script;
#                      used to create unique working directory names
#
# temp_dir: Java temp directory to use for jobs
#
#
# Author: David Roazen <droazen@broadinstitute.org>
#

if [ $# -ne 5 ]
then
    echo "Usage: $0 test_class_suffix job_queue timeout bamboo_build_number temp_dir"
    echo "Example test class suffixes: UnitTest, IntegrationTest, PipelineTest"
    exit 1
fi

TEST_CLASS_SUFFIX="$1"
JOB_QUEUE="$2"
GLOBAL_TIMEOUT="$3"
BUILD_NUMBER="$4"
TEMP_DIR="$5"

BAMBOO_CLONE=`pwd`
BAMBOO_BUILD_DIRECTORY=`basename "${BAMBOO_CLONE}"`
# name of the bamboo build directory + the build number can serve as a unique id for this run
BAMBOO_BUILD_ID="${BAMBOO_BUILD_DIRECTORY}-${BUILD_NUMBER}"
TEST_ROOT_WORKING_DIR="/humgen/gsa-scr1/gsa-engineering/parallel_tests_working_directory/${BAMBOO_BUILD_ID}"
TEST_CLONE="${TEST_ROOT_WORKING_DIR}/test_clone"
IVY_CACHE="${TEST_ROOT_WORKING_DIR}/ivy_cache"

JOB_OUTPUT_DIR="${TEST_ROOT_WORKING_DIR}/job_output"
JOB_MEMORY="4"
# We check job status every JOB_POLL_INTERVAL seconds
JOB_POLL_INTERVAL=30


# Kill any outstanding jobs
shutdown_jobs() {
    bjobs -p -r -P "${BAMBOO_BUILD_ID}" 2> /dev/null | grep -v JOBID | awk '{ print $1; }' | xargs bkill
}

# Print job output to stdout so that it gets recorded in bamboo's logs
echo_job_output() {
    cd "${JOB_OUTPUT_DIR}"
    for job_output in *.out
    do
        echo ""
        echo "-----------------------------------------------"
        echo "${job_output}:"
        echo "-----------------------------------------------"
        cat "${job_output}"
    done
}

# Setup working environment:

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

# Compile everything ONCE. All instances of the test suite will share this build:

cd "${TEST_CLONE}"
ant clean test.compile -Divy.home="${IVY_CACHE}" -Djava.io.tmpdir="${TEMP_DIR}"

if [ $? -ne 0 ]
then
    echo "$0: test.compile failed"
    exit 1
fi

# Submit one job per test class:

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

# Monitor jobs for completion. If they don't complete within the global timeout, terminate them:

ELAPSED_SECONDS=0
while [ "${ELAPSED_SECONDS}" -lt "${GLOBAL_TIMEOUT}" ]
do
    sleep "${JOB_POLL_INTERVAL}" 
    ELAPSED_SECONDS=`expr "${ELAPSED_SECONDS}" + "${JOB_POLL_INTERVAL}"`

    NUM_OUTSTANDING_JOBS=`bjobs -p -r -P "${BAMBOO_BUILD_ID}" 2> /dev/null | grep -v JOBID | wc -l`

    # If all jobs are done, save the test results back into bamboo's working directory,
    # and print the job output logs to stdout so that bamboo will save them in its log
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

        echo_job_output

        # TODO: delete job working directory before exiting (keep it for now for debugging purposes)
        exit 0
    fi

    echo "Outstanding jobs: ${NUM_OUTSTANDING_JOBS}"
done

echo "$0: Timeout of ${GLOBAL_TIMEOUT} seconds reached before jobs completed, giving up"
shutdown_jobs
exit 1


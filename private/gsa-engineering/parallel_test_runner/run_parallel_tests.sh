#!/bin/bash
#
# Simple implementation of class-level test suite parallelism. Dispatches one job
# to the farm per test class, monitors jobs for completion, and saves the aggregated
# results for display in Bamboo.
#
# Usage: run_parallel_tests.sh job_queue timeout bamboo_build_number temp_dir test_class_suffixes...
#
# Arguments:
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
# test_class_suffixes: run tests from test classes ending in these suffixes
#                      (eg., UnitTest, IntegrationTest, etc.)
#                      can specify an arbitrary number of suffixes
#
# Author: David Roazen <droazen@broadinstitute.org>
#

if [ $# -lt 5 ]
then
    echo "Usage: $0 job_queue timeout bamboo_build_number temp_dir test_class_suffixes..." 1>&2
    echo "Example test class suffixes: UnitTest, IntegrationTest, PipelineTest" 1>&2
    exit 1
fi

# Record the start time so that we can later record how long this run took
START_TIMESTAMP=`date '+%s'`

JOB_QUEUE="$1"
GLOBAL_TIMEOUT="$2"
BUILD_NUMBER="$3"
TEMP_DIR="$4"

# shift the args so that the variable number of trailing test_class_suffix args start at $1
shift 4

BAMBOO_CLONE=`pwd`
BAMBOO_BUILD_DIRECTORY=`basename "${BAMBOO_CLONE}"`
# name of the bamboo build directory + the build number can serve as a unique id for this run
BAMBOO_BUILD_ID="${BAMBOO_BUILD_DIRECTORY}-${BUILD_NUMBER}"
GLOBAL_PARALLEL_TESTS_DIR="/humgen/gsa-hpprojects/GATK/testing/parallel_tests_working_directories"
TEST_ARCHIVE_DIR="${GLOBAL_PARALLEL_TESTS_DIR}/archive"
TEST_ROOT_WORKING_DIR="${GLOBAL_PARALLEL_TESTS_DIR}/${BAMBOO_BUILD_ID}"
TEST_CLONE="${TEST_ROOT_WORKING_DIR}/test_clone"
IVY_CACHE="${TEST_ROOT_WORKING_DIR}/ivy_cache"
LOG_FILE="/humgen/gsa-hpprojects/GATK/testing/logs/parallel_tests_runtime.log"

JOB_RUNNER_DIR="${TEST_ROOT_WORKING_DIR}/job_runners"
JOB_OUTPUT_DIR="${TEST_ROOT_WORKING_DIR}/job_output"
JOB_MEMORY="4"
# We check job status every JOB_POLL_INTERVAL seconds
JOB_POLL_INTERVAL=30


# Kill any outstanding jobs
shutdown_jobs() {
    echo "$0: shutting down jobs for ${BAMBOO_BUILD_ID}"
    bjobs -A -J "${BAMBOO_BUILD_ID}" 2> /dev/null | grep -v JOBID | awk '{ print $1; }' | xargs bkill
    echo "$0: done shutting down jobs for ${BAMBOO_BUILD_ID}"
}

# Archive the working directory used by this parallel test suite run
# Deleting is too time-consuming -- old working directories can be deleted later by a cron job
archive_working_dir() {
    echo "$0: archiving test working directory"
    cd "${BAMBOO_CLONE}"
    mv "${TEST_ROOT_WORKING_DIR}" "${TEST_ARCHIVE_DIR}"

    if [ $? -eq 0 ]
    then
        echo "$0: done archiving test working directory"
    else
        echo "$0: failed to archive test working directory" 1>&2
    fi
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

# Record a log entry for this run including the Bamboo build ID, exit status (COMPLETED or TIMED_OUT),
# total runtime in seconds, and seconds spent waiting for farm jobs to complete
write_log_entry() {
    EXIT_STATUS="$1"
    END_TIMESTAMP=`date '+%s'`
    TOTAL_RUNTIME_IN_SECONDS=`expr "${END_TIMESTAMP}" - "${START_TIMESTAMP}"`

    printf "%s\t%s\t%d\t%d\n" "${BAMBOO_BUILD_ID}" "${EXIT_STATUS}" "${TOTAL_RUNTIME_IN_SECONDS}" "${SECONDS_WAITING_FOR_JOBS}" >> "${LOG_FILE}"
}


# Setup working environment:

if [ -d "${TEST_ROOT_WORKING_DIR}" ]
then
    rm -rf "${TEST_ROOT_WORKING_DIR}"
fi

mkdir "${TEST_ROOT_WORKING_DIR}"
mkdir "${IVY_CACHE}"
mkdir "${JOB_RUNNER_DIR}"
mkdir "${JOB_OUTPUT_DIR}"

echo "$0: Cloning ${BAMBOO_CLONE} into ${TEST_CLONE}"
git clone "${BAMBOO_CLONE}" "${TEST_CLONE}"

if [ $? -ne 0 ]
then
    echo "$0: failed to clone ${BAMBOO_CLONE} into ${TEST_CLONE}" 1>&2
    exit 1
fi

# Compile everything ONCE. All instances of the test suite will share this build:

cd "${TEST_CLONE}"
ant clean test.compile -Divy.home="${IVY_CACHE}" -Djava.io.tmpdir="${TEMP_DIR}"

if [ $? -ne 0 ]
then
    echo "$0: test.compile failed" 1>&2
    exit 1
fi

NUM_JOBS=0

# Create a separate shell script for each test class to be run

while [ $# -gt 0 ]
do
    TEST_CLASS_SUFFIX="$1"
    echo "Preparing job runners for ${TEST_CLASS_SUFFIX} test classes"
    NUM_CLASSES_THIS_SUFFIX=0

    for test_class_file in `find . -name "*${TEST_CLASS_SUFFIX}.class"`
    do
        test_class=`basename "${test_class_file}" ".class"`
        NUM_JOBS=`expr "${NUM_JOBS}" + 1`
        NUM_CLASSES_THIS_SUFFIX=`expr ${NUM_CLASSES_THIS_SUFFIX} + 1`

        echo "#!/bin/bash" > "${JOB_RUNNER_DIR}/${NUM_JOBS}.sh"
        echo "ant runtestonly \
              -Dsingle=${test_class} \
              -Divy.home=${IVY_CACHE} \
              -Djava.io.tmpdir=${TEMP_DIR} \
              > ${JOB_OUTPUT_DIR}/${test_class}.job.out 2>&1" \
              >> "${JOB_RUNNER_DIR}/${NUM_JOBS}.sh"

        chmod +x "${JOB_RUNNER_DIR}/${NUM_JOBS}.sh"
    done

    echo "Found ${NUM_CLASSES_THIS_SUFFIX} ${TEST_CLASS_SUFFIX} test classes"

    # move on to next test class suffix
    shift 1
done

# Submit all jobs at once as a single job array:

echo "Dispatching jobs for ${BAMBOO_BUILD_ID}"

bsub -P "${BAMBOO_BUILD_ID}" \
     -q "${JOB_QUEUE}" \
     -R "rusage[mem=${JOB_MEMORY}] select[tmp>100]" \
     -o "MASTER.job.out" \
     -J "${BAMBOO_BUILD_ID}[1-${NUM_JOBS}]" \
     "${JOB_RUNNER_DIR}/\${LSB_JOBINDEX}.sh"

if [ $? -ne 0 ]
then
    echo "Failed to dispatch job array" 1>&2
    exit 1
fi

echo "Dispatched ${NUM_JOBS} jobs to job queue ${JOB_QUEUE} as job array ${BAMBOO_BUILD_ID}[1-${NUM_JOBS}]"

# Monitor jobs for completion. If they don't complete within the global timeout, terminate them:

SECONDS_WAITING_FOR_JOBS=0
while [ "${SECONDS_WAITING_FOR_JOBS}" -lt "${GLOBAL_TIMEOUT}" ]
do
    sleep "${JOB_POLL_INTERVAL}" 
    SECONDS_WAITING_FOR_JOBS=`expr "${SECONDS_WAITING_FOR_JOBS}" + "${JOB_POLL_INTERVAL}"`

    NUM_OUTSTANDING_JOBS=`bjobs -J "${BAMBOO_BUILD_ID}" 2> /dev/null | grep -v JOBID | wc -l`

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
        archive_working_dir
        write_log_entry "COMPLETED"

        exit 0
    fi

    echo "--------------------------------------------------------------------------------"
    echo "Unfinished Jobs / Total Jobs: ${NUM_OUTSTANDING_JOBS}/${NUM_JOBS}"
    echo ""
    bjobs -A -J "${BAMBOO_BUILD_ID}"
    echo "--------------------------------------------------------------------------------"
done

echo "$0: Timeout of ${GLOBAL_TIMEOUT} seconds reached before jobs completed, giving up" 1>&2
write_log_entry "TIMED_OUT"

shutdown_jobs

# ok to spend time deleting the working dir in this case, since we've already spent an excessive
# amount of time on this run
cd "${BAMBOO_CLONE}" && rm -rf "${TEST_ROOT_WORKING_DIR}"

exit 1


#!/bin/bash
#
# Simple implementation of class-level test suite parallelism. Dispatches one job
# to the farm per test class, monitors jobs for completion, and saves the aggregated
# results for display in Bamboo.
#
# Usage: run_parallel_tests.sh job_queue timeout bamboo_build_number test_class_suffixes...
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
# test_class_suffixes: run tests from test classes ending in these suffixes
#                      (eg., UnitTest, IntegrationTest, etc.)
#                      can specify an arbitrary number of suffixes
#
# Author: David Roazen <droazen@broadinstitute.org>
#

if [ $# -lt 4 ]
then
    echo "Usage: $0 job_queue timeout bamboo_build_number test_class_suffixes..." 1>&2
    echo "Example test class suffixes: UnitTest, IntegrationTest, PipelineTest" 1>&2
    exit 1
fi

# Record the start time so that we can later record how long this run took
START_TIMESTAMP=`date '+%s'`

JOB_QUEUE="$1"
GLOBAL_TIMEOUT="$2"
BUILD_NUMBER="$3"

# shift the args so that the variable number of trailing test_class_suffix args start at $1
shift 3

BAMBOO_CLONE=`pwd`
BAMBOO_BUILD_DIRECTORY=`basename "${BAMBOO_CLONE}"`
# name of the bamboo build directory + the build number can serve as a unique id for this run
BAMBOO_BUILD_ID="${BAMBOO_BUILD_DIRECTORY}-${BUILD_NUMBER}"
COMPILE_IVY_CACHE="${BAMBOO_CLONE}/ivy_cache"
COMPILE_TEMP_DIR="${BAMBOO_CLONE}/tmp"
GLOBAL_PARALLEL_TESTS_DIR="/broad/hptmp/gsabamboo/parallel_tests_working_directories"
TEST_ARCHIVE_DIR="${GLOBAL_PARALLEL_TESTS_DIR}/archive"
TEST_ROOT_WORKING_DIR="${GLOBAL_PARALLEL_TESTS_DIR}/${BAMBOO_BUILD_ID}"
TEST_CLONE="${TEST_ROOT_WORKING_DIR}/test_clone"
TEST_TEMP_DIR="${TEST_ROOT_WORKING_DIR}/tmp"
LOG_FILE="/local/gsa-engineering/parallel_tests_logs/parallel_tests_runtime.log"

JOB_RUNNER_DIR="${TEST_ROOT_WORKING_DIR}/job_runners"
JOB_OUTPUT_DIR="${TEST_ROOT_WORKING_DIR}/job_output"
JOB_MEMORY="4"
# We check job status every JOB_POLL_INTERVAL seconds
JOB_POLL_INTERVAL=30

# We attempt to cd into each of these directories before each job in an effort to avoid automount failures
declare -a AUTOMOUNT_DIR_LIST=( "/seq/references/" \
                                "/humgen/1kg/reference/" \
                                "/humgen/gsa-hpprojects/" \
                              )
# Seconds to wait after attempting to trigger automount
AUTOMOUNT_TRIGGER_DELAY=10


# Print the names of any test classes that have not yet finished running
print_unfinished_test_class_names() {
    UNFINISHED_TEST_CLASSES="Unfinished Test Classes: "

    for UNFINISHED_TEST_CLASS_ID in `bjobs -J "${BAMBOO_BUILD_ID}" | grep -v "DONE" | grep -v "EXIT" | grep -v "JOBID" | sed 's/^.*\[//g' | sed 's/\].*$//g'`
    do
        UNFINISHED_TEST_CLASS=`cat "${JOB_RUNNER_DIR}/${UNFINISHED_TEST_CLASS_ID}.sh" | grep runtestonly | awk '{ print $3; }' | awk -F'=' '{ print $2; }'`
        UNFINISHED_TEST_CLASSES="${UNFINISHED_TEST_CLASSES} ${UNFINISHED_TEST_CLASS}"
    done

    echo ""
    echo "${UNFINISHED_TEST_CLASSES}"
}

# Print the names of any test classes that exited with an error
print_exited_test_class_names() {
    EXITED_TEST_CLASSES="ERROR In Test Classes: "

    for EXITED_TEST_CLASS_ID in `bjobs -a -J "${BAMBOO_BUILD_ID}" | grep -v "JOBID" | grep "EXIT" | sed 's/^.*\[//g' | sed 's/\].*$//g'`
    do
        EXITED_TEST_CLASS=`cat "${JOB_RUNNER_DIR}/${EXITED_TEST_CLASS_ID}.sh" | grep runtestonly | awk '{ print $3; }' | awk -F'=' '{ print $2; }'`
        EXITED_TEST_CLASSES="${EXITED_TEST_CLASSES} ${EXITED_TEST_CLASS}"
    done

    echo ""
    echo "${EXITED_TEST_CLASSES}"
}

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

# Detect the case where there were no actual test failures, but one or more jobs exited with an error
check_for_non_test_related_job_failures() {
    if ! grep 'BUILD FAILED' "${JOB_OUTPUT_DIR}"/*
    then
        if [ "${JOBS_EXITED_WITH_ERROR}" -ne 0 ]
        then
            echo "No test failures, but ${JOBS_EXITED_WITH_ERROR} jobs exited with an error. Probably a temporary farm glitch." 1>&2

            archive_working_dir
            write_log_entry "FARM_GLITCH"
            exit 1
        fi
    fi
}

create_automount_triggers_command_line() {
    AUTOMOUNT_COMMAND_LINE=""

    for automount_dir in ${AUTOMOUNT_DIR_LIST[@]}
    do
        AUTOMOUNT_COMMAND_LINE="${AUTOMOUNT_COMMAND_LINE} cd \"${automount_dir}\";"
    done

    AUTOMOUNT_COMMAND_LINE="${AUTOMOUNT_COMMAND_LINE} sleep ${AUTOMOUNT_TRIGGER_DELAY}"
    echo "${AUTOMOUNT_COMMAND_LINE}"
}

# Setup working environment:

if [ ! -d "${GLOBAL_PARALLEL_TESTS_DIR}" ]
then
    mkdir -p "${GLOBAL_PARALLEL_TESTS_DIR}"
fi

if [ ! -d "${TEST_ARCHIVE_DIR}" ]
then
    mkdir -p "${TEST_ARCHIVE_DIR}"
fi

if [ -d "${TEST_ROOT_WORKING_DIR}" ]
then
    rm -rf "${TEST_ROOT_WORKING_DIR}"
fi

mkdir "${COMPILE_IVY_CACHE}"
mkdir "${COMPILE_TEMP_DIR}"
mkdir "${TEST_ROOT_WORKING_DIR}"
mkdir "${TEST_TEMP_DIR}"
mkdir "${JOB_RUNNER_DIR}"
mkdir "${JOB_OUTPUT_DIR}"

# Introduce a delay to make sure bamboo clone has been fully written to disk
sleep 5

# Compile everything ONCE. All instances of the test suite will share this build:
ant clean test.compile -Divy.home="${COMPILE_IVY_CACHE}" -Djava.io.tmpdir="${COMPILE_TEMP_DIR}"

if [ $? -ne 0 ]
then
    echo "$0: test.compile failed" 1>&2
    exit 1
fi

echo "$0: Copying ${BAMBOO_CLONE} into ${TEST_CLONE}"
cp -r "${BAMBOO_CLONE}" "${TEST_CLONE}"

if [ $? -ne 0 ]
then
    echo "$0: failed to copy ${BAMBOO_CLONE} into ${TEST_CLONE}" 1>&2
    exit 1
fi

cd "${TEST_CLONE}"

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
        echo "`create_automount_triggers_command_line`" >> "${JOB_RUNNER_DIR}/${NUM_JOBS}.sh"
        echo "cd ${TEST_CLONE}" >> "${JOB_RUNNER_DIR}/${NUM_JOBS}.sh"
        echo "ant runtestonly \
              -Dsingle=${test_class} \
              -Djava.io.tmpdir=${TEST_TEMP_DIR} \
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

    JOBS_EXITED_WITH_ERROR=`bjobs -a -A -J "${BAMBOO_BUILD_ID}" 2> /dev/null | grep -v JOBID | awk '{ print $8; }'`
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
        check_for_non_test_related_job_failures
        archive_working_dir
        write_log_entry "COMPLETED"

        exit 0
    fi

    echo "--------------------------------------------------------------------------------"
    echo "Unfinished Jobs / Total Jobs: ${NUM_OUTSTANDING_JOBS}/${NUM_JOBS}"
    echo ""
    bjobs -A -J "${BAMBOO_BUILD_ID}"

    if [ "${JOBS_EXITED_WITH_ERROR}" -gt 0 ]
    then
        print_exited_test_class_names
    fi

    if [ "${NUM_OUTSTANDING_JOBS}" -lt 10 ]
    then
        print_unfinished_test_class_names
    fi

    echo "--------------------------------------------------------------------------------"
done

echo "$0: Timeout of ${GLOBAL_TIMEOUT} seconds reached before jobs completed, giving up" 1>&2
write_log_entry "TIMED_OUT"

shutdown_jobs

# ok to spend time deleting the working dir in this case, since we've already spent an excessive
# amount of time on this run
cd "${BAMBOO_CLONE}" && rm -rf "${TEST_ROOT_WORKING_DIR}"

exit 1


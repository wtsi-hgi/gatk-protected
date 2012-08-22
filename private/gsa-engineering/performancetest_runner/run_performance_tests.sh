#!/bin/bash

TEST_NAME="$1"
IVY_CACHE_DIR="$2"
TEMP_DIR="$3"
BAMBOO_CLONE=`pwd`
PERFORMANCE_TEST_ROOT_DIR="/humgen/gsa-scr1/gsa-engineering/performance_test_working_directories"
PERFORMANCE_TEST_WORKING_DIR="${PERFORMANCE_TEST_ROOT_DIR}/${TEST_NAME}"
PERFORMANCE_TEST_QSCRIPT="private/scala/qscript/org/broadinstitute/sting/queue/qscripts/performance/GATKPerformanceOverTime.scala"
RESOURCES_DIR="/home/unix/depristo/work/oneOffProjects/gatkPerformanceOverTime/resources"
JOB_QUEUE="gsa-engineering"
ANALYSIS_SCRIPT="private/R/GATKPerformanceOverTime.R"
WEB_DIR="/humgen/gsa-hpprojects/performance"
ANNOUNCE_MAILING_LIST="gsamembers@broadinstitute.org"
EMAIL_COMMAND="/usr/sbin/sendmail -t"

if [ $# -ne 3 ]
then
    echo "$0: Usage: $0 test_name ivy_cache_dir temp_dir"
    exit 1
fi

if [ ! -d "${PERFORMANCE_TEST_ROOT_DIR}" ]
then
    echo "$0: ${PERFORMANCE_TEST_ROOT_DIR} does not exist. Please create it as user gsa-engineering."
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

if [ ! -d "${WEB_DIR}" ]
then
    echo "$0: Web directory ${WEB_DIR} for publishing performance reports does not exist. Please create it as user gsa-engineering."
    exit 1
fi

if [ ! -d "${RESOURCES_DIR}" ]
then
    echo "$0: Resources directory ${RESOURCES_DIR} does not exist. Please create and populate it."
    exit 1
fi

if [ -d "${PERFORMANCE_TEST_WORKING_DIR}" ]
then
    rm -rf "${PERFORMANCE_TEST_WORKING_DIR}"
fi

echo "$0: Cloning ${BAMBOO_CLONE} into ${PERFORMANCE_TEST_WORKING_DIR}"
git clone "${BAMBOO_CLONE}" "${PERFORMANCE_TEST_WORKING_DIR}"

if [ $? -ne 0 ]
then
    echo "$0: Failed to clone ${BAMBOO_CLONE} into ${PERFORMANCE_TEST_WORKING_DIR}"
    exit 1
fi

cd "${PERFORMANCE_TEST_WORKING_DIR}"

ant clean gsalib -Djava.io.tmpdir=${TEMP_DIR} -Divy.home=${IVY_CACHE_DIR}

if [ $? -ne 0 ]
then
    echo "$0: Failed to build gsalib."
    exit 1
fi

ant -Djava.io.tmpdir=${TEMP_DIR} -Divy.home=${IVY_CACHE_DIR} 

if [ $? -ne 0 ]
then
    echo "$0: Failed to build the GATK/Queue."
    exit 1
fi

java -Djava.io.tmpdir=${TEMP_DIR} \
-jar dist/Queue.jar \
-S "${PERFORMANCE_TEST_QSCRIPT}" \
-myJarFile dist/GenomeAnalysisTK.jar \
-resources "${RESOURCES_DIR}" \
-iterations 3 \
-qsub \
-jobQueue "${JOB_QUEUE}" \
-run \
-startFromScratch \
-l DEBUG

if [ $? -ne 0 ]
then
    echo "$0: Error running GATKPerformanceOverTime QScript." 
    exit 1
fi

# this gets the latest report file generated
REPORT=`ls -ltr *.jobreport.txt | tail -n 1 | awk '{print $9}'`
REPORT_BASENAME=`basename "${REPORT}" .txt`
echo "$0: Using report file ${REPORT}"

# debug
cp "${REPORT}" "${PERFORMANCE_TEST_ROOT_DIR}/debug"
echo "$0: copied ${REPORT} into debug location ${PERFORMANCE_TEST_ROOT_DIR}/debug"

# run the analysis script
# export R_LIBS_USER=`pwd`/public/R
Rscript "${ANALYSIS_SCRIPT}" "${REPORT}" "${REPORT}.performance_over_time.pdf"

if [ $? -ne 0 ]
then
    echo "$0: Error running GATKPerformanceOverTime R script." 
    exit 1
fi

# copy reports to web-accessible location, indexed by date and time
TIMESTAMP=`date '+%Y-%m-%d_%H:%M:%S'`
mkdir ${WEB_DIR}/${TIMESTAMP} && chmod 755 ${WEB_DIR}/${TIMESTAMP} && cp ${REPORT_BASENAME}* ${WEB_DIR}/${TIMESTAMP} && chmod 644 ${WEB_DIR}/${TIMESTAMP}/${REPORT_BASENAME}*

if [ $? -ne 0 ]
then
    echo "$0: Failed to copy ${REPORT_BASENAME}* files into web directory ${WEB_DIR}/${TIMESTAMP}" 
    exit 1
fi

cat <<-EOF | ${EMAIL_COMMAND}
	To: ${ANNOUNCE_MAILING_LIST}
	Subject: GATK performance log for ${TIMESTAMP} posted 
	Content-type: text/plain
	
	The GATK performance log for ${TIMESTAMP} has been posted:
	
	http://iwww.broadinstitute.org/gsa/performance/${TIMESTAMP}
	
	EOF

if [ $? -ne 0 ]
then
	echo "$0: Failed to send email to ${ANNOUNCE_MAILING_LIST} announcing the posting of the new performance logs."
	exit 1
fi

echo "$0: performance test suite run was successful"
exit 0


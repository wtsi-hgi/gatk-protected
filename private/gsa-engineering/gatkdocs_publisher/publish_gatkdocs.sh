#!/bin/bash
#
# Generate a new set of gatkdocs and upload them to a live location on the gsaweb server
#

WEB_SERVER="gsaweb"
GATKDOCS_LIVE_DIR="/local/htdocs/gatk/documentation/tooldocs"
GATKDOCS_STAGING_DIR="/local/htdocs/staging/gatkdocs"
GATKDOCS_LOCAL_DIR="target/gatkdocs"
TEMP_MAVEN_REPO="tmp_mvn_repo"
TEMP_DIR="tmp"

mkdir "${TEMP_DIR}"

mvn clean && \
mvn install "-Dmaven.repo.local=${TEMP_MAVEN_REPO}" '-P!private,!queue' -Ddisable.queue "-Djava.io.tmpdir=${TEMP_DIR}" -Dgatkdocs.extension=php && \
mvn site "-Dmaven.repo.local=${TEMP_MAVEN_REPO}" '-P!private,!queue' -Ddisable.queue "-Djava.io.tmpdir=${TEMP_DIR}" -Dgatkdocs.extension=php

if [ $? -ne 0 ]
then
    echo "Failed to generate gatkdocs"
    exit 1
fi

echo "Uploading gatkdocs to ${WEB_SERVER}"

ssh "${WEB_SERVER}" "rm -rf ${GATKDOCS_STAGING_DIR} && mkdir ${GATKDOCS_STAGING_DIR}" && \
rsync -rvtz "${GATKDOCS_LOCAL_DIR}"/* "${WEB_SERVER}:${GATKDOCS_STAGING_DIR}" && \
ssh "${WEB_SERVER}" "chmod -R 775 ${GATKDOCS_STAGING_DIR} && mv ${GATKDOCS_LIVE_DIR} ${GATKDOCS_LIVE_DIR}_old && mv ${GATKDOCS_STAGING_DIR} ${GATKDOCS_LIVE_DIR} && rm -rf ${GATKDOCS_LIVE_DIR}_old"

if [ $? -ne 0 ]
then
    echo "Failed to upload gatkdocs to web server"
    exit 1
fi


exit 0

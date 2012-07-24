#!/bin/bash

IWWW_DIR="gsa-stage:/local/software/apache2/htdocs/gatk/gatkdocs_private"
WWW_DIR="gsa-stage:/local/software/apache2/htdocs/gatk/gatkdocs"
STAGING_DIR="/local/gsa-engineering/gatkdocs_publisher/staging_area"
BAMBOO_CLONE=`pwd`

if [ $# -ne 1 ]
then
    echo "Usage: $0 <repository_name>"
    exit 1
fi

REPOSITORY="$1"

if [ "${REPOSITORY}" == "unstable" -o "${REPOSITORY}" == "stable" ]
then
    DESTINATION_DIR="${IWWW_DIR}/${REPOSITORY}" 
    GATKDOCS_INCLUDE_HIDDEN="-Dgatkdocs.include.hidden=true"
elif [ "${REPOSITORY}" == "release" -o "${REPOSITORY}" == "literelease" ]
then
    DESTINATION_DIR="${WWW_DIR}/${REPOSITORY}"
    GATKDOCS_INCLUDE_HIDDEN=""
else
    echo "$0: Invalid repository specified: ${REPOSITORY}"
    exit 1
fi

STAGING_CLONE="${STAGING_DIR}/${REPOSITORY}"

if [ -d "${STAGING_CLONE}" ]
then
    rm -rf "${STAGING_CLONE}"
fi

git clone "${BAMBOO_CLONE}" "${STAGING_CLONE}"

if [ $? -ne 0 ]
then
    echo "$0: Failed to clone ${BAMBOO_CLONE} into ${STAGING_CLONE}"
    exit 1
fi

IVY_CACHE_DIR="${STAGING_CLONE}/ivy_cache"
mkdir "${IVY_CACHE_DIR}"

cd "${STAGING_CLONE}"

printf "$0: working dir is: %s\n" `pwd`
printf "$0: destination dir is: %s\n" "${DESTINATION_DIR}"

if [ "${REPOSITORY}" == "release" ]
then
    rm -rf private
elif [ "${REPOSITORY}" == "literelease" ]
then
    rm -rf protected private
fi

ant clean && ant gatkdocs ${GATKDOCS_INCLUDE_HIDDEN} -Divy.home=${IVY_CACHE_DIR}

if [ $? -ne 0 ]
then
    echo "$0: Failed to generate gatkdocs"
    exit 1
fi

if [ ! -d "${DESTINATION_DIR}" ]
then
    mkdir "${DESTINATION_DIR}"
fi

if [ -d "${DESTINATION_DIR}_new" ] 
then
    rm -rf "${DESTINATION_DIR}_new"
fi

if [ -d "${DESTINATION_DIR}_old" ]
then
    rm -rf "${DESTINATION_DIR}_old"
fi

# Try to minimize the window of time in which the docs are unavailable
# as we replace the old docs with the new docs:

cp -r gatkdocs "${DESTINATION_DIR}_new" && \
rsync "${DESTINATION_DIR}" "${DESTINATION_DIR}_old" && \
rsync "${DESTINATION_DIR}_new" "${DESTINATION_DIR}" && \
rm -rf "${DESTINATION_DIR}_old"

if [ $? -ne 0 ]
then
    echo "Failed to copy gatkdocs into ${DESTINATION_DIR}"
    exit 1
fi

exit 0


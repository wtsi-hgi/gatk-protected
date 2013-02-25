#!/bin/bash

IWWW_DIR="gsa-web:/local/software/apache2/htdocs/gatkdocs_private"
WWW_DIR="gsa-web:/local/software/apache2/htdocs/gatk/gatkdocs"
STAGING_DIR="/local/gsa-engineering/gatkdocs_publisher/staging_area"
BAMBOO_CLONE=`pwd`

if [ $# -ne 1 ]
then
    echo "Usage: $0 <repository_name>"
    exit 1
fi

REPOSITORY="$1"

if [ "${REPOSITORY}" == "unstable" ] 
then 
	DESTINATION_DIR="${IWWW_DIR}"
    GATKDOCS_INCLUDE_HIDDEN="-Dgatkdocs.include.hidden=true"
elif [ "${REPOSITORY}" == "release" ]
then
	DESTINATION_DIR="${WWW_DIR}"
    GATKDOCS_INCLUDE_HIDDEN=""
else
	echo "No need to update gatkdocs for ${REPOSITORY} -- everything is fine"
	exit 0
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
fi

ant clean && ant gatkdocs ${GATKDOCS_INCLUDE_HIDDEN} -Divy.home=${IVY_CACHE_DIR}

if [ $? -ne 0 ]
then
    echo "$0: Failed to generate gatkdocs"
    exit 1
fi

rsync -rvtz --delete gatkdocs/* ${DESTINATION_DIR} 

if [ $? -ne 0 ]
then
    echo "Failed to copy gatkdocs into ${DESTINATION_DIR}"
    exit 1
fi

exit 0


#!/bin/bash

GATK_PUBLIC_REPO_URL="git@github.com:broadgsa/gatk.git"
GATK_PROTECTED_REPO_URL="git@github.com:broadgsa/gatk-protected.git"

if [ $# -ne 1 ]
then
    echo "Usage: $0 release_target"
    exit 1
fi

if [ "$1" == "public" ]
then
    TARGET_REPO="${GATK_PUBLIC_REPO_URL}"
elif [ "$1" == "protected" ]
then
    TARGET_REPO="${GATK_PROTECTED_REPO_URL}"
else
    echo "$0: Invalid release target: $1"
    exit 1
fi

git remote add github "${TARGET_REPO}"

NEW_HEAD=`git rev-parse HEAD`
REMOTE_HEAD=`git ls-remote --heads github | grep refs/heads/master | awk '{print $1}'`

if [ "${NEW_HEAD}" == "${REMOTE_HEAD}" ]
then
    echo "$0: Release repository ${TARGET_REPO} up-to-date with respect to Stable, no source release performed (HEAD = ${REMOTE_HEAD})"
    exit 0
fi

git push --tags github master:master 

if [ $? -ne 0 ]
then
    echo "$0: SOURCE RELEASE FAILED: Failed to push into ${RELEASE_URL}"
    exit 1
fi

echo "$0: SOURCE RELEASE SUCCEEDED: Successfully updated Release repository from latest Stable (new HEAD = ${NEW_HEAD})"

exit 0


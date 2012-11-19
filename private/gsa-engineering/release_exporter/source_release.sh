#!/bin/bash

RELEASE_URL="git@github.com:broadgsa/gatk.git"

git remote add github "${RELEASE_URL}"

NEW_HEAD=`git rev-parse HEAD`
REMOTE_HEAD=`git ls-remote --heads github | grep refs/heads/master | awk '{print $1}'`

if [ "${NEW_HEAD}" == "${REMOTE_HEAD}" ]
then
    echo "$0: Release repository up-to-date with respect to Stable, no source release performed (HEAD = ${REMOTE_HEAD})"
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


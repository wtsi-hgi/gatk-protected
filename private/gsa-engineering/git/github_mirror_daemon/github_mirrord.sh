#!/bin/bash
#
# Daemon process that keeps a directory of github mirrors up-to-date.
# Updates the mirrors every POLLING_INTERVAL seconds.
#
# Author: David Roazen
#

MIRROR_DIR="/local/gsa-engineering/github_mirrors" 
POLLING_INTERVAL=30

if [ ! -d "${MIRROR_DIR}" ]
then
    echo "$0: ${MIRROR_DIR} does not exist."
    exit 1
fi

cd "${MIRROR_DIR}"

while true
do
    for mirror in `find . -mindepth 1 -maxdepth 1 -type d`
    do
        cd "${mirror}"
        git remote update origin
        git remote prune origin
        cd ..
    done

    sleep "${POLLING_INTERVAL}"
done


#!/bin/bash
# 
# Respawn the github mirror daemon if necessary. 
#
# Author: David Roazen
#

DAEMON="/local/gsa-engineering/github_mirrors/github_mirrord.sh"
DAEMON_NAME="github_mirrord"
DAEMON_USER="gsa-engineering"

if ! ps -U "${DAEMON_USER}" -o comm | grep "${DAEMON_NAME}"
then
    echo "Github mirror daemon not found, restarting ${DAEMON}"
    nohup "${DAEMON}" > /dev/null &
    sleep 5 
else
    echo "Github mirror daemon running, no action required"
fi

exit 0


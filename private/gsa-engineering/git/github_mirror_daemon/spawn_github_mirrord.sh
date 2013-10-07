#!/bin/bash
# 
# Respawn the github mirror daemon if necessary. 
#
# Author: David Roazen
#

DAEMON="/local/gsa-engineering/github_mirrors/github_mirrord.sh"
DAEMON_NAME="github_mirrord"
DAEMON_USER="gsa-engineering"

if ! ps -U "${DAEMON_USER}" -o comm | grep -q "${DAEMON_NAME}"
then
    nohup "${DAEMON}" > /dev/null &
    sleep 5 
fi

exit 0


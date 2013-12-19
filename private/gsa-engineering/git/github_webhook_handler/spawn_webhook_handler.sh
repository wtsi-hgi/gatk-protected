#!/bin/bash
# 
# Respawn the github webhook handler if necessary. 
#
# Author: David Roazen
#

HANDLER="/humgen/gsa-scr1/gsa-engineering/git/github_webhook_handler/handle_webhooks.sh"
HANDLER_NAME="handle_webhooks"
HANDLER_USER="gsa-engineering"

if ! ps -U "${HANDLER_USER}" -o comm | grep -q "${HANDLER_NAME}"
then
    nohup "${HANDLER}" > /dev/null &
    sleep 5 
fi

exit 0


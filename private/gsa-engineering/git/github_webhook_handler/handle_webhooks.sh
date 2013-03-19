#!/bin/bash
#
# Handles logged github webhook events by feeding them to the post-receive email script.
# Run as a cron job at frequent intervals.
#
# Author: David Roazen
#

LOG_DIR="/humgen/gsa-scr1/gsa-engineering/git/github_webhook_logs"
MIRROR_DIR="/humgen/gsa-scr1/gsa-engineering/git/github_webhook_handler/mirrors"
EMAIL_SCRIPT="/humgen/gsa-scr1/gsa-engineering/git/email_scripts/post-receive-email"

if [ ! -d "${LOG_DIR}" ]
then
    echo "$0: ${LOG_DIR} does not exist"
    exit 1
fi

if [ ! -d "${MIRROR_DIR}" ]
then
    echo "$0: ${MIRROR_DIR} does not exist"
    exit 1
fi

if [ ! -x "${EMAIL_SCRIPT}" ]
then
    echo "$0: ${EMAIL_SCRIPT} does not exist or is not executable"
    exit 1
fi

handle_webhook_activity() {
    shopt -s nullglob

    for mirror in ${MIRROR_DIR}/*
    do
        mirror_basename=`basename $mirror`
        mirror_log_dir="${LOG_DIR}/${mirror_basename}"

        echo "Checking repo ${mirror_basename}"

        test -d "${mirror}" || continue
        test -d "${mirror_log_dir}" || continue

        cd "${mirror}"

        echo "Now in ${mirror}"
        echo "Checking logs in ${mirror_log_dir}"

        for log in ${mirror_log_dir}/*
        do
            echo "Found log $log"

            git remote update origin

            if [ $? -eq 0 ]
            then
                entry=`head -1 "$log"`
                echo "$entry" | ${EMAIL_SCRIPT}

                if [ $? -eq 0 ]
                then
                    rm -f "$log"
                else
                    echo "Error sending email for log entry $log"
                fi
            else
                echo "Error updating local mirror ${mirror}"
            fi
        done
    done
}

handle_webhook_activity
exit 0

#!/bin/bash
# 
# Respawn the Bamboo server if necessary. 
#
# Author: David Roazen
#

BAMBOO_ROOT="/local/gsa-engineering/bamboo"
BAMBOO_HOME="/local/gsa-engineering/bamboo/atlassian-bamboo-4.4.1"
BAMBOO_PID_FILE="${BAMBOO_HOME}/bamboo.pid"
BAMBOO_STARTUP_SCRIPT="${BAMBOO_ROOT}/startup_bamboo_with_correct_env.sh"
BAMBOO_SERVER_CLASS="com.atlassian.bamboo.server.Server"
EMERGENCY_MAIL_ADDRESS="gsa-engineering@broadinstitute.org"
SENDMAIL="/usr/sbin/sendmail -t"

BAMBOO_PID=`cat "${BAMBOO_PID_FILE}" | head -1 | awk '{ print $1; }'`

send_mail() {
    MESSAGE="$1"

    cat <<-EOF | ${SENDMAIL}
	From: gsa-engineering@broadinstitute.org
	To: ${EMERGENCY_MAIL_ADDRESS}
	Subject: $MESSAGE 
	Content-type: text/plain

	$0: $MESSAGE 
	EOF
}

if ! ps -p "${BAMBOO_PID}" -f | grep -q "${BAMBOO_SERVER_CLASS}"
then
    cd "${BAMBOO_ROOT}"
    bash "${BAMBOO_STARTUP_SCRIPT}" 
    
    if [ $? -ne 0 ]
    then
        send_mail "BAMBOO SERVER DOWN ON GSA4, UNABLE TO RESTART"
        exit 1
    fi

    send_mail "RESTARTED BAMBOO SERVER ON GSA4"
fi

exit 0


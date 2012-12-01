#!/bin/bash

ANNOUNCE_MAILING_LIST="gsa-announce@broadinstitute.org"
EMAIL_COMMAND="/usr/sbin/sendmail -t"

GATK_FULL_VERSION=`git describe --long`
if [ $? -ne 0 ]
then
	echo "$0: Failed to determine new GATK version number (git describe failed)."
	exit 1
fi
GATK_VERSION=`echo "${GATK_FULL_VERSION}" | awk -F '-' '{ print $1 "-" $2; }'`

cat <<-EOF | ${EMAIL_COMMAND}
	From: Broad GSA <gsa-engineering@broadinstitute.org>
	To: ${ANNOUNCE_MAILING_LIST}
	Subject: GATK version ${GATK_VERSION} released
	Content-type: text/plain
	
	Version ${GATK_VERSION} of the GATK has just been released.
	
	To download it, please visit http://www.broadinstitute.org/gatk/download 
	
	
	If this is a major (eg., 1.0 -> 2.0) or minor (eg., 2.0 -> 2.1) update,
	release notes will be posted on our website:
	
	http://gatkforums.broadinstitute.org/discussions/tagged/release-notes
	
	If this is a bug fix to an existing release (eg., 2.1-1 -> 2.1-2), you 
	can see what has changed by examining the commit logs in our github
	repository:
	
	https://github.com/broadgsa/gatk/commits/master
	
	
	For documentation and tutorials, please visit the main GATK website:
	
	http://www.broadinstitute.org/gatk/
	
	Or, if our documentation doesn't answer your question, please post your question on our support forum:
	
	http://gatkforums.broadinstitute.org/
	
	
	Regards,
	
	The GSA Team
	
	EOF

if [ $? -ne 0 ]
then
	echo "$0: Failed to send email to ${ANNOUNCE_MAILING_LIST} announcing the release of GATK version ${GATK_VERSION}"
	exit 1
fi

echo "$0: Successfully sent email to ${ANNOUNCE_MAILING_LIST} announcing the release of GATK version ${GATK_VERSION}"
exit 0



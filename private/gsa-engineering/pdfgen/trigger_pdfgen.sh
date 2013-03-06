#!/bin/bash

PDFGEN_URL="http://www.broadinstitute.org/gatk/private/pdfgen?token=j1a9b7a9u0m4a1n2n&version="
GIT_VERSION=`git describe --long | awk -F'-' '{ print $1 "-" $2; }'`

wget -q -O /dev/null "${PDFGEN_URL}${GIT_VERSION}"

if [ $? -ne 0 ]
then
    echo "Failed to trigger the pdfgen script"
    exit 1
fi

exit 0


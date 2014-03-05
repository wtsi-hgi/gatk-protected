#!/bin/bash

WORKING_DIR=`pwd`
PRIVATE_TESTDATA="${WORKING_DIR}/private/gatk-private/src/test/resources"

declare -a LINK_DIR_LIST=( "private/gatk-private" \
                           "protected/gatk-protected" \
                           "public/gatk-framework" \
                           "public/gatk-package" \
                           "public/queue-framework" \
                           "public/queue-package" \
                           "${WORKING_DIR}" \
                         )

for dir in ${LINK_DIR_LIST[@]}
do
    cd "${WORKING_DIR}"
    cd "${dir}"

    mkdir private
    ln -fs "${PRIVATE_TESTDATA}" private/testdata
done

exit 0


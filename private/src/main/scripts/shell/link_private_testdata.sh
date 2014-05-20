#!/bin/bash

WORKING_DIR=`pwd`
PRIVATE_TESTDATA="${WORKING_DIR}/private/gatk-tools-private/src/test/resources"

declare -a LINK_DIR_LIST=( "private/gatk-package-internal" \
                           "private/gatk-queue-extensions-internal" \
                           "private/gatk-queue-package-internal" \
                           "private/gatk-tools-private" \
                           "protected/gatk-package-distribution" \
                           "protected/gatk-queue-extensions-distribution" \
                           "protected/gatk-queue-package-distribution" \
                           "protected/gatk-tools-protected" \
                           "public/gatk-engine" \
                           "public/gatk-queue" \
                           "public/gatk-queue-extensions-public" \
                           "public/gatk-tools-public" \
                           "public/gatk-utils" \
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


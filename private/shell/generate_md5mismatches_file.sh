#!/bin/bash
#
# Given a Bamboo integration test log file as an argument, creates an
# md5mismatches.txt file in the current directory suitable for feeding
# to the private/python/updateMD5s.py script for a bulk MD5 update
#

if [ $# -ne 1 ]
then
    echo "Usage: $0 bamboo_integration_test_log_file"
    exit 1
fi

BAMBOO_LOG="$1"

printf "expected\tobserved\ttest\n" > md5mismatches.txt
grep "MD5 mismatch" "${BAMBOO_LOG}" | awk '{ print $10 "\t" $8 "\t" "foo" }' | sort | uniq >> md5mismatches.txt
exit 0


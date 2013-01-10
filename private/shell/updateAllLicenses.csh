#!/bin/tcsh

# This script will update all files in the GATK with the license information
# in the files living in the license directory. This script is meant to be
# run manually (and hopefully only once). To run it, you must be in the
# $GATK/private/shell directory. If you want to run it from a different directory
# you will have to update the relative paths on this file and aldo the ones on
# updateLicenseForFile.csh
#
# author: Mauricio Carneiro
# date: 1/9/13

echo "";
echo "Updating all files in the GATK with the new licenses";
echo '- Run with "-debug" flag to list all the files being updated';
echo "";

echo "Processing private directory..."
find ../../private/java   -name "*.java" -exec updateLicenseForFile.csh {} ../../licensing/private_license.txt $1 \;   ;

echo "Processing public directory..."
find ../../public/java    -name "*.java" -exec updateLicenseForFile.csh {} ../../licensing/public_license.txt $1 \;    ;

echo "Processing protected directory..."
find ../../protected/java -name "*.java" -exec updateLicenseForFile.csh {} ../../licensing/protected_license.txt $1 \; ;


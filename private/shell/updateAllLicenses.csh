#!/bin/tcsh

# This script will update all files in the GATK with the license information
# in the files living in the license directory. This script is meant to be
# run manually (and hopefully only once). To run it, you must be in the
# $GATK root directory. If you want to run it from a different directory
# you will have to update the relative paths on this file and aldo the ones on
# updateLicenseForFile.csh
#
# author: Mauricio Carneiro
# date: 1/9/13

echo "";
echo "Updating all files in the GATK with the new licenses";
echo '- Run with "-debug" flag to list all the files being updated';
echo "";

echo "Processing java private directory..."
find private   -name "*.java" -exec private/shell/updateLicenseForFile.csh {} licensing/private_license.txt $1 \;   ;

echo "Processing java public directory..."
find public    -name "*.java" -exec private/shell/updateLicenseForFile.csh {} licensing/public_license.txt $1 \;    ;

echo "Processing java protected directory..."
find protected -name "*.java" -exec private/shell/updateLicenseForFile.csh {} licensing/protected_license.txt $1 \; ;



echo "Processing scala private directory..."
find private    -name "*.scala" -exec private/shell/updateLicenseForFile.csh {} licensing/private_license.txt $1 \;   ;

echo "Processing scala public directory..."
find public     -name "*.scala" -exec private/shell/updateLicenseForFile.csh {} licensing/public_license.txt $1 \;    ;

echo "Processing scala protected directory..."
find protected  -name "*.scala" -exec private/shell/updateLicenseForFile.csh {} licensing/protected_license.txt $1 \; ;


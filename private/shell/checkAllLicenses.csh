#!/bin/tcsh

# This script will CHECK all files in the GATK for the license information
# it is meant to be run by bamboo as a sanity check that our repo contains the
# correct license information in all files.
#
# script must be run from the $GIT_ROOT
#
# author: Mauricio Carneiro
# date: 1/9/13

echo "";
echo "Checking all licenses in the GATK.. ";

find private   -name "*.java"  |  python private/python/CheckLicense.py licensing/private_license.txt   
if ($status != "0") then
	exit(1)
endif

find public    -name "*.java"  |  python private/python/CheckLicense.py licensing/public_license.txt    
if ($status != "0") then
	exit(1)
endif

find protected -name "*.java"  |  python private/python/CheckLicense.py licensing/protected_license.txt 
if ($status != "0") then
	exit(1)
endif

find private   -name "*.scala" |  python private/python/CheckLicense.py licensing/private_license.txt   
if ($status != "0") then
	exit(1)
endif

find public    -name "*.scala" |  python private/python/CheckLicense.py licensing/public_license.txt    
if ($status != "0") then
	exit(1)
endif

find protected -name "*.scala" |  python private/python/CheckLicense.py licensing/protected_license.txt;


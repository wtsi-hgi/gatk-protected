#!/bin/tcsh

if ("$2" == "-debug") then
    echo $1
endif

if (!($1:e == "java" || $1:e == "scala")) then
    echo "Only java or scala files should be updated for license. This is a $1:e file"
    exit 1
endif

set isPrivate = `echo "$1" | sed -n '/private/p'`;
set isProtected = `echo "$1" | sed -n '/protected/p'`;
set isPublic = `echo "$1" | sed -n '/public/p'`;

if ($isPrivate != "") then
    set licenseFile = licensing/private_license.txt
else if ($isProtected != "") then
    set licenseFile = licensing/protected_license.txt
else if ($isPublic != "") then
    set licenseFile = licensing/public_license.txt
else
    echo "is not in public, private or protected and is a source file. Please place it in the appropriate directory so we can choose the appropriate license";
    exit 1;
endif

private/python/ParseLicense.py $1 $licenseFile > $1.newlicense && mv $1.newlicense $1

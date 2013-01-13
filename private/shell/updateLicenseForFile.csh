#!/bin/tcsh

if ("$3" == "-debug") then
    echo $1
endif

cat $1 | ../python/ParseLicense.py $2 > $1.newlicense && mv $1.newlicense $1;

#!/bin/tcsh
set ret = `grep "package" $1`

if ( "$ret" == "" ) then
	echo "$1"
endif



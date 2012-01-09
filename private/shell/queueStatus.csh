#!/bin/tcsh

echo "\nGSA queue usage"
bjobs -u all -q gsa | awk '$2 !~ "USER" {print $2}' | sort | uniq -c

echo "\nGeneral computing resources"
bqueues gsa week hour

echo "\nFH jobs"
bjobs -u gsaadm

echo "\nFile system status"
df -h /humgen/* /broad/shptmp 

# what's the status of all of the gsa hosts
echo "\nGSA host status"
bhosts gsahosts


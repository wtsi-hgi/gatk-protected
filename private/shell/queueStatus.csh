#!/bin/tcsh

echo "\nGSA queue usage"
bjobs -u all -q gsa | awk '$2 !~ "USER" {print $2}' | sort | uniq -c

echo "\nGSAFolk usage"
busers -w gsafolk

echo "\nGSA folk"
bugroup -l gsafolk

echo "\nGeneral computing resources"
bqueues gsa week hour

echo "\nFH jobs"
bjobs -u gsaadm | wc -l

echo "\nFile system status"
df -h /humgen/* /broad/shptmp 

# what's the status of all of the gsa hosts
echo "\nGSA host status"
bhosts gsahosts


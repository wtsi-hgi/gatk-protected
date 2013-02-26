#!/bin/tcsh

set java = $1
shift
echo $java -cp build/java/classes:build/java/testclasses:lib/testng-5.14.1.jar:lib/log4j-1.2.16.jar:lib/commons-lang-2.5.jar org.broadinstitute.sting.utils.nanoScheduler.NanoSchedulerUnitTest $*

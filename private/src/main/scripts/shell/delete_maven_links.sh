#!/bin/bash
#
# Delete symlinks to testdata/qscripts created by maven. This is necessary
# in cases where "mvn clean" fails to delete these links itself and exits
# with an error.
#
# Should be run from the root directory of your git clone.
#

find . -type l -name testdata -exec rm '{}' ';'
find . -type l -name qscript -exec rm '{}' ';'

exit 0


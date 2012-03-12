#!/bin/tcsh

setenv CLEAN_DIR ~/dev/cleanSourceCopies

foreach type (unstable stable gatk.git)
pushd $CLEAN_DIR/$type > /dev/null
git pull -q > /dev/null
setenv X '$0'
echo "type $type"
git rev-list master | xargs git describe 
popd > /dev/null
end


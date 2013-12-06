#!/bin/tcsh

setenv CLEAN_DIR /local/gsa-engineering/cron_clones

foreach type (unstable stable gatk.git)
pushd $CLEAN_DIR/$type > /dev/null
git pull -q > /dev/null
setenv X '$0'
echo "type $type"
git rev-list master | xargs git describe --long
popd > /dev/null
end


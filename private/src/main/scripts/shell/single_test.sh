#!/bin/sh

# Run a single passed in test on jars in an already installed local repository.
# Primarily for use by build systems such as bamboo to parallelize tests.

test_name=$1
local_repo=${2-"local_repo"}

test_property="it.test"
test_disabled="test"
test_type="commit"

if [[ "${test_name}" =~ .*UnitTest ]] ; then
    # UnitTests run with -Dtest=MyUnitTest instead of it.test
    test_property="test"
    test_disabled="it.test"
elif [[ "${test_name}" =~ .*PipelineTest ]] ; then
    # Pipeline tests are historically defined as separate from "commit"
    test_type="pipeline"
fi

echo mvn verify \
  -Dmaven.repo.local=${local_repo} \
  -Dmaven.javadoc.skip=true \
  -Dsting.generate-gatk-extensions.skipped=true \
  -Dsting.jar.phase=none \
  -Dsting.unpack.phase=none \
  -Dsting.shade.phase=none \
  -Dsting.package${test_type}tests.skipped=false \
  -D${test_disabled}=disabled \
  -D${test_property}=${test_name}

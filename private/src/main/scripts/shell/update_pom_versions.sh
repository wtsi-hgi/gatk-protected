#!/bin/sh

version=${1-`git describe --abbrev=0`-SNAPSHOT}


echo "updating root" && \
cd public/sting-root && \
mvn versions:set -DnewVersion=${version} -DgenerateBackupPoms=false && \

# Gave up on versions:set for invoked-only pom non-child pom, due to perplexing errors:
# a) "Non-resolvable parent POM: Could not find artifact"
# b) "Project version is inherited from parent"

echo "updating package tests" && \
cd ../package-tests && \
sed -i -e '/^        <artifactId>sting-root<\/artifactId>$/ {
    N
    s@\(        <artifactId>sting-root</artifactId>\n        <version>\).*\(</version>\)@\1'${version}'\2@
}' pom.xml && \

echo "updating external example" && \
cd ../external-example && \
sed -i -e 's@<sting.version>.*\</sting.version>@<sting.version>'${version}'</sting.version>@' pom.xml && \

echo "done"

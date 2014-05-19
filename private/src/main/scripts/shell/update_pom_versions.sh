#!/bin/sh

if [ "$1" == "" ]; then
    echo "Usage: $0 <version>" >&2
    exit 1
fi

version=${1-`git describe --abbrev=0`-SNAPSHOT}

echo "updating root" && \
cd public/gatk-root && \
mvn versions:set -DnewVersion=${version} -DgenerateBackupPoms=false && \

# Gave up on versions:set for invoked-only pom non-child pom, due to perplexing errors:
# a) "Non-resolvable parent POM: Could not find artifact"
# b) "Project version is inherited from parent"

cd .. && \
echo "updating package tests and Vector Pair HMM" && \
sed -i -e '/^        <artifactId>gatk-root<\/artifactId>$/ {
    N
    s@\(        <artifactId>gatk-root</artifactId>\n        <version>\).*\(</version>\)@\1'${version}'\2@
}' package-tests/pom.xml VectorPairHMM/pom.xml && \

echo "updating external example" && \
sed -i -e 's@<gatk.version>.*\</gatk.version>@<gatk.version>'${version}'</gatk.version>@' external-example/pom.xml && \

echo "done"

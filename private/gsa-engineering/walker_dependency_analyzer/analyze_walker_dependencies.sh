#!/bin/bash
#
# Frontend script to the DependencyAnalyzer that takes a list of walker classes and a pair of
# git commit IDs, and determines which walkers have compile-time dependencies on classes changed
# between the two commits. Considers all classes in the */java/src directories, as well as classes
# within jar files that have changed.
#
# Outputs a Java properties file, with each key/value pair indicating whether a walker has dependencies
# on any of the changed classes. Example output file contents:
#
#     org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyper=true
#     org.broadinstitute.sting.gatk.walkers.haplotypecaller.HaplotypeCaller=true
#     org.broadinstitute.sting.gatk.walkers.readutils.PrintReads=false
#
# Prerequisites:
#     -Must be run from the root directory of a clone of one of our git repositories.
#     -The bcel and ant-apache-bcel jars must be installed in ~/.ant/lib
#
# Usage: private/gsa-engineering/walker_dependency_analyzer/analyze_walker_dependencies.sh walker_file output_file old_git_commit_id new_git_commit_id
#
#        walker_file: A file containing a list of the walker classes to analyze, one per line. Walker class names must
#                     be fully-qualified.
#                         Example:
#                         org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyper
#                         org.broadinstitute.sting.gatk.walkers.haplotypecaller.HaplotypeCaller
#                         org.broadinstitute.sting.gatk.walkers.readutils.PrintReads
#
#        output_file: File to which to write the results. Results will be formatted as discussed above.
#
#        old_git_commit_id:
#        new_git_commit_id: Valid git commit IDs from which to produce the diff of changed files.
#                           All java source files changed between old_git_commit_id and new_git_commit_id
#                           will be included in the analysis.
#
# Author: David Roazen
#


if [ $# -ne 4 ]
then
    echo "Usage: $0 walker_file output_file old_git_commit_id new_git_commit_id"
    exit 1
fi

WALKER_FILE="$1"
OUTPUT_FILE="$2"
OLD_GIT_COMMIT_ID="$3"
NEW_GIT_COMMIT_ID="$4"
DEPENDENCY_ANALYZER_LOCATION="private/gsa-engineering/walker_dependency_analyzer"
DEPENDENCY_ANALYZER_NAME="DependencyAnalyzer.xml"
CHANGED_CLASSES_FILE="analyze_walker_dependencies_changed_classes"

if [ ! -f "${WALKER_FILE}" ]
then
    echo "$0: walker file ${WALKER_FILE} does not exist"
    exit 1
fi

if [ -f "${OUTPUT_FILE}" ]
then
    echo "$0: output file ${OUTPUT_FILE} already exists, refusing to overwrite."
    exit 1
fi

git rev-parse --verify "${OLD_GIT_COMMIT_ID}" > /dev/null && git rev-parse --verify "${NEW_GIT_COMMIT_ID}" > /dev/null
if [ $? -ne 0 ]
then
    echo "$0: either ${OLD_GIT_COMMIT_ID} or ${NEW_GIT_COMMIT_ID} was an invalid git commit ID for this repository"
    exit 1
fi

if [ -f "${CHANGED_CLASSES_FILE}" ]
then
    echo "$0: Found stale ${CHANGED_CLASSES_FILE} file, removing"
    rm "${CHANGED_CLASSES_FILE}"
fi

# Include all changed classes in the */java/src directories
git diff --name-only "${OLD_GIT_COMMIT_ID}" "${NEW_GIT_COMMIT_ID}" | grep -E '^(public|protected|private)/java/src/' | cut -d'/' -f4- | sed 's/\.java$/\.class/g' >> "${CHANGED_CLASSES_FILE}"

# Also include all classes within jars that have changed
for changed_jar in `git diff --name-only "${OLD_GIT_COMMIT_ID}" "${NEW_GIT_COMMIT_ID}" | grep -E '\.jar$'`
do
    # Both the old and new versions of the jar will appear in the diff if the jar name has changed.
    # Make sure to only use the current jar (ie., the one that still exists):
    if [ -f "${changed_jar}" ]
    then
        jar -tvf "${changed_jar}" | awk '{ print $8; }' | grep -E '\.class$' >> "${CHANGED_CLASSES_FILE}"
    fi
done

ant clean package.gatk.all

if [ $? -ne 0 ]
then
    echo "$0: failed to compile and package the GATK in preparation for dependency analysis"
    exit 1
fi

cp "${DEPENDENCY_ANALYZER_LOCATION}/${DEPENDENCY_ANALYZER_NAME}" .

for walker_class in `cat "${WALKER_FILE}"`
do
    ant -f "${DEPENDENCY_ANALYZER_NAME}" \
        -Dwalker="${walker_class}" \
        -Dchanged.classes.file="${CHANGED_CLASSES_FILE}" \
        -Dresult.properties.file="${OUTPUT_FILE}"
done

echo "$0: dependency analysis complete, results are in ${OUTPUT_FILE}."
echo "Echoing results:"
cat "${OUTPUT_FILE}"

exit 0

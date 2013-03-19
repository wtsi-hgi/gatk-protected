#!/bin/bash
#
# Removes blacklisted files and directories from the git history of the clone it's run in. Must only be 
# run by the gsa-engineering user, and only on a clone of Stable.
#
# Author: David Roazen
#

# Directories to filter out of the git history for the specified release. IMPORTANT: If we ever change these lists
# (EXCEPT to blacklist a brand-new directory that didn't previously exist in the repository) we'll need to 
# rebuild the github repository from scratch, and everyone who has cloned it will need to delete their old 
# clones and re-clone it. This is due to the way history re-writing in git works: the commit ids all change 
# to reflect the new hashes of the directory contents and the new hashes of their ancestor commits, so if 
# you change the directories to be deleted you end up creating an entirely new git history unrelated to the 
# former history, making it impossible to update the repository through a simple "git push" and breaking
# synchronization for everyone. In other words, THINK CAREFULLY BEFORE CHANGING THIS VARIABLE!
PUBLIC_REPO_BLACKLIST="private protected"
PROTECTED_REPO_BLACKLIST="private"

GATK_PUBLIC_REPO_URL="git@github.com:broadgsa/gatk.git"
GATK_PROTECTED_REPO_URL="git@github.com:broadgsa/gatk-protected.git"

if [ $# -ne 1 ]
then
    echo "Usage: $0 release_target"
    exit 1
fi

if [ "$1" == "public" ]
then
    BLACKLIST="${PUBLIC_REPO_BLACKLIST}"
    GITHUB_URL="${GATK_PUBLIC_REPO_URL}"
elif [ "$1" == "protected" ]
then
    BLACKLIST="${PROTECTED_REPO_BLACKLIST}"
    GITHUB_URL="${GATK_PROTECTED_REPO_URL}"
else
    echo "$0: invalid release target $1"
    exit 1
fi

echo "$0: Prior to filtering, local HEAD = `git rev-parse HEAD`"

echo "$0: Filtering history in clone `pwd` to remove blacklisted directories"
git filter-branch --index-filter "git rm -rf --cached --ignore-unmatch ${BLACKLIST}" --tag-name-filter cat | grep -v '^rm ' 

if [ $? -ne 0 ]
then
    echo "$0: git filter-branch exited with an error"
    exit 1
fi

# Introduce a delay here in an attempt to reduce transient git errors
sleep 10

echo "$0: After filtering, local HEAD = `git rev-parse HEAD`"

CURRENT_GITHUB_HEAD=`git ls-remote --heads ${GITHUB_URL} | grep refs/heads/master | awk '{print $1}'`
git show ${CURRENT_GITHUB_HEAD} > /dev/null

if [ $? -ne 0 ]
then
    echo "$0: Unable to find current github HEAD (${CURRENT_GITHUB_HEAD}) within the rewritten local history, refusing to proceed!"
    exit 1
fi

echo "$0: Located the current github HEAD (${CURRENT_GITHUB_HEAD}) within the rewritten local history. Filtering succeeded!"

exit 0


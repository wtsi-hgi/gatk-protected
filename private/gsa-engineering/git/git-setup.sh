#!/bin/bash

if [ $# -ne 2 ]
then
    echo "Usage: $0 \"Full Name\" \"Email Address\""
    exit 1
fi

if [ -f ~/.gitconfig ] && grep -q '^\[url' ~/.gitconfig
then 
    cat ~/.gitconfig | grep -v '^\[url' | grep -v "insteadOf" > ~/.gitconfig_cleaned
    mv -f ~/.gitconfig_cleaned ~/.gitconfig
fi

USER_NAME="$1"
USER_EMAIL="$2"

STABLE_URL="git@github.com:broadinstitute/gsa-stable.git"
UNSTABLE_URL="git@github.com:broadinstitute/gsa-unstable.git"
PUBLIC_GITHUB_URL="git://github.com/broadgsa/"

git config --global user.name "${USER_NAME}"
git config --global user.email "${USER_EMAIL}"

git config --global "url.${STABLE_URL}.insteadOf" gsa:stable
git config --global "url.${UNSTABLE_URL}.insteadOf" gsa:unstable
git config --global "url.${PUBLIC_GITHUB_URL}.insteadOf" github:

git config --global remote.stable.url "${STABLE_URL}"
git config --global remote.stable.fetch '+refs/heads/*:refs/remotes/stable/*'

git config --global remote.unstable.url "${UNSTABLE_URL}"
git config --global remote.unstable.fetch '+refs/heads/*:refs/remotes/unstable/*'

git config --global push.default upstream

git config --global alias.fix-unstable '!git fetch unstable && git checkout -b unstable unstable/master && git merge -m "Merged bug fix from Stable into Unstable" master'

git config --global alias.push-unstable-fix '!git push unstable unstable:master && git checkout master && git branch -D unstable'

git config --global alias.diff-origin '!git fetch origin && git diff origin/master HEAD'
git config --global alias.diffstat-origin '!git fetch origin && git diff --stat=200,200 origin/master HEAD'
git config --global alias.web 'instaweb --httpd=webrick'

echo "Git is now set up and ready to use!"
exit 0


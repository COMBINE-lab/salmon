#!/bin/bash

if [ $# -eq 0 ]
then
    echo "No input arguments provided.  Usage is merge_into_develop.sh <branch_to_merge_in>"
    exit 1
fi

feature=$1

# from https://stackoverflow.com/questions/173919/is-there-a-theirs-version-of-git-merge-s-ours
# in case branchA is not our current branch
git checkout develop

# make merge commit but without conflicts!!
# the contents of 'ours' will be discarded later
git merge -s ours ${feature}

# make temporary branch to merged commit
git branch branchTEMP

# get contents of working tree and index to the one of branchB
git reset --hard ${feature}

# reset to our merged commit but 
# keep contents of working tree and index
git reset --soft branchTEMP

# change the contents of the merged commit
# with the contents of branchB
git commit --amend

# get rid off our temporary branch
git branch -D branchTEMP

# verify that the merge commit contains only contents of branchB
git diff HEAD ${feature}

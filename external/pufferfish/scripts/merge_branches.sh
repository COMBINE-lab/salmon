#!/bin/bash

target=""
source=""

PARAMS=""
while (( "$#" )); do
    case "$1" in
        -t|--merge-into)
            target=$2
            shift 2
            ;;
        -f|--merge-from)
            source=$2
            shift 2
            ;;
        --) # end argument parsing
            shift
            break
            ;;
        -*|--*=) # unsupported flags
        echo "Error: Unsupported flag $1" >&2
        exit 1
        ;;
        *) # preserve positional arguments
            PARAM="$PARAMS $1"
            shift
            ;;
    esac
done
# set positional arguments in their proper place
eval set -- "$PARAMS"

if [[ -z "${source// }" ]]; then
   echo "Source branch is empty" >&2
   exit 1
fi

if [[ -z "${target// }" ]]; then
   echo "Target branch is empty" >&2
   exit 1
fi

echo "merging ${source} into ${target} using branchTEMP"

# in case branchA is not our current branch
git checkout $target

# make merge commit but without conflicts!!
# the contents of 'ours' will be discarded later
git merge -s ours $source

# make temporary branch to merged commit
git branch branchTEMP

# get contents of working tree and index to the one of branchB
git reset --hard $source

# reset to our merged commit but 
# keep contents of working tree and index
git reset --soft branchTEMP

# change the contents of the merged commit
# with the contents of branchB
git commit --amend 

# get rid off our temporary branch
git branch -D branchTEMP

# verify that the merge commit contains only contents of branchB
git diff HEAD ${source}

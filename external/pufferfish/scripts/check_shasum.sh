#!/bin/bash

exists()
{
    command -v "$1" >/dev/null 2>&1
}

sum=$1
fname=$2

hashcheck=""
if exists sha256sum; then
	hashcheck="sha256sum"
elif exists shasum; then
	hashcheck="shasum -a256"
else
	unset hashcheck
fi

if [ -z "${hashcheck-}" ]; then
    echo "Couldn't find shasum command; can't verify contents of ${fname}";
else
    echo "${sum}  ${fname}" | ${hashcheck} -c - || { echo "${fname} did not match expected SHA256! Exiting."; exit 1; }
fi



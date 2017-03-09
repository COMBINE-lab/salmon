#!/bin/bash

exists()
{
    command -v "$1" >/dev/null 2>&1
}

sum=$1
fname=$2

if exists sha256sum; then
    echo "${sum}  ${fname}" | sha256sum -c - || { echo "${fname} did not match expected SHA256! Exiting."; exit 1; }
else
    echo "Couldn't find shasum command; can't verify contents of ${fname}";
fi



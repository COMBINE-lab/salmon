#!/bin/bash
GTF=$1
DIR="$( dirname -- "$1"  )"
OUT="$DIR/gid_tid.txt"
#DIR="$(dirname -- "$(realpath -- "$1")")"
awk -F "\t" '$3 == "transcript" {print $9}' $GTF | tr -d ";\"" | awk '{print $2"\t"$4}' > $OUT
echo "gid_tid.txt written in "$DIR

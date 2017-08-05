#!/bin/bash
FASTA=$1
DIR="$( dirname -- "$1"  )"
OUT="$DIR/transcript_clean.fasta"
#DIR="$(dirname -- "$(realpath -- "$1")")"
sed -e '/^>/ s/|.*//' $FASTA > $OUT
#awk -F "\t" '$3 == "transcript" {print $9}' $GTF | tr -d ";\"" | awk '{print $2"\t"$4}' > $OUT
echo "clean fasta file written here: "$OUT

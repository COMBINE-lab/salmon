#!/bin/bash

cmd="$@"
file1=`echo $cmd | grep -Po -- '-1\s+\K[^ ]+'`
new_cmd=`echo $cmd | sed 's/-1[[:space:]][[:graph:]]\+//'`
file2=`echo $cmd | grep -Po -- '-2\s+\K[^ ]+'`
new_cmd=`echo $new_cmd | sed 's/-2[[:space:]][[:graph:]]\+//'`

tmpdir=$(mktemp -d)
echo "TEMPDIR is $tmpdir"
# Cleanup on exit
trap 'rm -rf "$tmpdir"' EXIT INT TERM HUP
p1="$tmpdir/p1.fa"
p2="$tmpdir/p2.fa"
mkfifo $p1
mkfifo $p2

./wrapper $file1 $file2 >> $p1 2>> $p2 &

echo "Running command [${new_cmd} -1 $p1 -2 $p2]"
${new_cmd} -1 $p1 -2 $p2

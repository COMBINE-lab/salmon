#!/bin/bash

cmd="$@"
base=`echo $cmd | grep -Po -- '-b\s+\K[^ ]+'`
new_cmd=`echo $cmd | sed 's/-b[[:space:]][[:graph:]]\+//'`

tmpdir=$(mktemp -d)
echo "TEMPDIR is $tmpdir"
# Cleanup on exit
trap 'rm -rf "$tmpdir"' EXIT INT TERM HUP
p1="$tmpdir/p1.fa"
p2="$tmpdir/p2.fa"
mkfifo $p1
mkfifo $p2

i1=`ls $base*I1*`
wrapper <(cat $base*I1*) <(cat $base*RA*) >> $p1 2>> $p2 &

echo "Running command [${new_cmd} -1 $p1 -2 $p2 -r $i1]"
${new_cmd} -1 $p1 -2 $p2 -r $i1

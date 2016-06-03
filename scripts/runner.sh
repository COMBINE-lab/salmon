#!/bin/bash

cmd="$@"
interleaved_file=`echo $cmd | sed -n 's/.*--interleaved\s\+\(\S\+\)\s\+.*/\1/p'`

if [ -z "$interleaved_file" ]
then
    #Run normally in this branch
    ${@}
else
   new_cmd=`echo $cmd | sed 's/--interleaved\s\+\S\+\s\+//'`
   tmpdir=$(mktemp -d)
   # Cleanup on exit
   trap 'rm -rf "$tmpdir"' EXIT INT TERM HUP
   p1="$tmpdir/p1.fq"
   p2="$tmpdir/p2.fq"
    mkfifo $p1
    mkfifo $p2
    # The following interleaved to split conversion is courtesy of
    # https://gist.github.com/nathanhaigh/3521724
    (paste - - - - - - - - | tee >(cut -f 1-4 | tr '\t' '\n' > $p1) | cut -f 5-8 | tr '\t' '\n' > $p2) < $interleaved_file &
    echo "Running command [${new_cmd} -1 $p1 -2 $p2]"
    ${new_cmd} -1 $p1 -2 $p2
fi
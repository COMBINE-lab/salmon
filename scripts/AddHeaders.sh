#!/bin/bash

usage()
{
cat << EOF
usage: $0 options

This script adds a header to all files matching the provided pattern in the given directory

OPTIONS:
   -h      Show this message
   -l      license file 
   -d      directory to search
   -p      file pattern to match
EOF
}

license=
pattern=
directory=
while getopts "hl:p:d:" OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         l)
             license=$OPTARG
             ;;
         p)
             pattern=$OPTARG
             ;;
         d)
             directory=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

echo "Prepending ${license} to files with pattern ${pattern} in directory ${directory}"

for i in ${directory}/${pattern}
do
  if ! grep -q Copyright $i
  then
    cat ${license} $i >$i.new && mv $i.new $i
  fi
done

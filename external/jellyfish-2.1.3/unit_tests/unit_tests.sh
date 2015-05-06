#! /bin/sh

cd tests
. ./compat.sh

if [ -n "$VALGRIND" ]; then
    exec valgrind ${DIR}/test_all
else
    exec ${DIR}/test_all
fi

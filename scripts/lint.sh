#!/bin/bash
set -e

# cmakelint
# https://github.com/richq/cmake-lint

if ! which cmakelint; then
    echo "Install cmakelint." 1>&2
    exit 1
fi
cmake_files=$(find . -name CMakeLists.txt -o -name "*.cmake")
cmakelint --config=.cmakelintrc ${cmake_files}

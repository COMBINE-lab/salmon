#!/bin/bash
set -e

source /hbb_shlib/activate

set -x

export CFLAGS="-g -O2 -I/hbb_shlib/include"
export CXXFLAGS="-g -O2 -I/hbb_shlib/include"

CPATH=`pwd`
echo "[Drone build] current path : ${CPATH}"
echo "[Drone build] making build directory"
mkdir build
cd build

echo "[Drone build] cmake configuration"
cmake ..

echo "[Drone build] making pufferfish and installing locally (this could take a while)"
make

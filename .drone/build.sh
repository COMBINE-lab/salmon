#!/bin/bash
source /hbb_exe/activate

set -e

CPATH=`pwd`
echo "[Drone build] current path : ${CPATH}"
echo "[Drone build] making build directory"

mkdir build
cd build

echo "[Drone build] cmake configuration"

cmake -DDO_QUIET_MAKE=TRUE -DBOOST_ROOT=/usr ..

echo "[Drone build] making salmon and installing locally (this could take a while)"

make -j8 -s install

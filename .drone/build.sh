#!/bin/bash
source /hbb_shlib/activate

set -e

CPATH=`pwd`
echo "[Drone build] current path : ${CPATH}"
echo "[Drone build] making build directory"

mkdir build
cd build

echo "[Drone build] cmake configuration"

cmake -DDO_QUIET_MAKE=TRUE -DBOOST_ROOT=/usr -DNO_IPO=TRUE -DCMAKE_BUILD_TYPE=RELEASE ..

echo "[Drone build] making salmon and installing locally (this could take a while)"

make -j8 -s install

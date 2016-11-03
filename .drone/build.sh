#!/bin/bash
set -e

CPATH=`pwd`
echo "[Drone build] current path : ${CPATH}"
echo "[Drone build] making build directory"

mkdir build
cd build

echo "[Drone build] cmake configuration"

cmake ..

echo "[Drone build] making salmon and installing locally (this could take a while)"

make install

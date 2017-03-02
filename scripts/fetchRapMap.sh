#!/bin/bash

exists()
{
  command -v "$1" >/dev/null 2>&1
}

CURR_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
EXTERNAL_DIR=${CURR_DIR}/../external
INSTALL_DIR=${CURR_DIR}/../external/install

if [ -d ${EXTERNAL_DIR}/RapMap ] ; then
    rm -fr ${EXTERNAL_DIR}/RapMap
fi

if [ -d ${INSTALL_DIR}/include/rapmap ] ; then
    rm -fr ${INSTALL_DIR}/include/rapmap
fi

if [ -d ${INSTALL_DIR}/src/rapmap ] ; then
    rm -fr ${INSTALL_DIR}/src/rapmap
fi

mkdir -p ${EXTERNAL_DIR}
#curl -k -L https://github.com/COMBINE-lab/RapMap/archive/salmon-v0.8.0.zip -o ${EXTERNAL_DIR}/rapmap.zip
curl -k -L https://github.com/COMBINE-lab/RapMap/archive/develop-salmon.zip -o ${EXTERNAL_DIR}/rapmap.zip

if exists shasum; then
  echo "d2462af6f66a4ee95a92add65b0663d21b507b2c  ${EXTERNAL_DIR}/rapmap.zip" | shasum -a1 -c - || { echo "rapmap.zip did not match expected SHA1! Exiting."; exit 1; }
else
  echo "Couldn't find shasum command; can't verify contents of downloaded RapMap";
fi


rm -fr ${EXTERNAL_DIR}/RapMap
unzip ${EXTERNAL_DIR}/rapmap.zip -d ${EXTERNAL_DIR}
mv ${EXTERNAL_DIR}/RapMap-develop-salmon ${EXTERNAL_DIR}/RapMap
#mv ${EXTERNAL_DIR}/RapMap-salmon-v0.8.0 ${EXTERNAL_DIR}/RapMap

mkdir -p ${INSTALL_DIR}/include/rapmap
mkdir -p ${INSTALL_DIR}/src/rapmap

rm ${EXTERNAL_DIR}/RapMap/src/xxhash.c
rm ${EXTERNAL_DIR}/RapMap/include/xxhash.h

cp -r ${EXTERNAL_DIR}/RapMap/external/libdivsufsort.zip ${EXTERNAL_DIR}
cp -r ${EXTERNAL_DIR}/RapMap/src/*.c ${INSTALL_DIR}/src/rapmap
cp -r ${EXTERNAL_DIR}/RapMap/src/*.cpp ${INSTALL_DIR}/src/rapmap
cp -r ${EXTERNAL_DIR}/RapMap/include/tclap ${INSTALL_DIR}/include/rapmap
cp -r ${EXTERNAL_DIR}/RapMap/include/*.h ${INSTALL_DIR}/include/rapmap
cp -r ${EXTERNAL_DIR}/RapMap/include/*.hpp ${INSTALL_DIR}/include/rapmap

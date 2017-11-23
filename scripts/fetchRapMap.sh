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
#curl -k -L https://github.com/COMBINE-lab/RapMap/archive/salmon-v0.8.2.zip -o ${EXTERNAL_DIR}/rapmap.zip
curl -k -L https://github.com/COMBINE-lab/RapMap/archive/salmon-v0.9.0.zip -o ${EXTERNAL_DIR}/rapmap.zip
#curl -k -L https://github.com/COMBINE-lab/RapMap/archive/develop-salmon.zip -o ${EXTERNAL_DIR}/rapmap.zip

hashcheck=""
if exists sha256sum; then
	hashcheck="sha256sum"
elif exists shasum; then
	hashcheck="shasum -a256"
else
	unset hashcheck
fi

if [ -z "${hashcheck-}" ]; then
    echo "Couldn't find shasum command; can't verify contents of downloaded RapMap";
else
    #echo "1691f4bca2b604f05f36772ae45faf0842ab4809843df770bd10366a5cfd6822  ${EXTERNAL_DIR}/rapmap.zip" | ${hashcheck} -c - || { echo "rapmap.zip did not match expected SHA1! Exiting."; exit 1; }
    #echo "5e0fa601f97c7b5f3bfa3216755d20a957d3f5e761f1128566a1441145cb60d7  ${EXTERNAL_DIR}/rapmap.zip" | ${hashcheck} -c - || { echo "rapmap.zip did not match expected SHA1! Exiting."; exit 1; }
    echo "e3dd29e3f79ad93f655c2262f428688d3dd81462d5ba3b7191dc1b16cbe5489f  ${EXTERNAL_DIR}/rapmap.zip" | ${hashcheck} -c - || { echo "rapmap.zip did not match expected SHA1! Exiting."; exit 1; }
    #echo "not testing sha in develop branch"
fi


rm -fr ${EXTERNAL_DIR}/RapMap
unzip ${EXTERNAL_DIR}/rapmap.zip -d ${EXTERNAL_DIR}
#mv ${EXTERNAL_DIR}/RapMap-develop-salmon ${EXTERNAL_DIR}/RapMap
mv ${EXTERNAL_DIR}/RapMap-salmon-v0.9.0 ${EXTERNAL_DIR}/RapMap

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
cp -r ${EXTERNAL_DIR}/RapMap/include/sparsepp ${INSTALL_DIR}/include/rapmap

#!/bin/bash
set -eu -o pipefail

exists()
{
  command -v "$1" >/dev/null 2>&1
}

CURR_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
EXTERNAL_DIR=${CURR_DIR}/../external
INSTALL_DIR=${CURR_DIR}/../external/install

if [ -d ${EXTERNAL_DIR}/pufferfish ] ; then
    rm -fr ${EXTERNAL_DIR}/pufferfish
fi

if [ -d ${INSTALL_DIR}/include/pufferfish ] ; then
    rm -fr ${INSTALL_DIR}/include/pufferfish
fi

if [ -d ${INSTALL_DIR}/src/pufferfish ] ; then
    rm -fr ${INSTALL_DIR}/src/pufferfish
fi

SVER=salmon-v1.4.0
#SVER=develop
#SVER=sketch-mode

EXPECTED_SHA256=059207e8d3134060ed70595e53f4189954c9e5edfaa6361b46304f55d1b71bc7

mkdir -p ${EXTERNAL_DIR}
curl -k -L https://github.com/COMBINE-lab/pufferfish/archive/${SVER}.zip -o ${EXTERNAL_DIR}/pufferfish.zip

hashcheck=""
if exists sha256sum; then
	hashcheck="sha256sum"
elif exists shasum; then
	hashcheck="shasum -a256"
else
	unset hashcheck
fi

if [ -z "${hashcheck-}" ]; then
    echo "Couldn't find shasum command; can't verify contents of downloaded pufferfish";
else
    if [[ $SVER != develop ]]; then
        echo "${EXPECTED_SHA256}  ${EXTERNAL_DIR}/pufferfish.zip" | ${hashcheck} -c - || { echo "pufferfish.zip did not match expected SHA1! Exiting."; exit 1; }
    else
        echo "not testing sha since pulling from develop"
    fi
fi


rm -fr ${EXTERNAL_DIR}/pufferfish
unzip ${EXTERNAL_DIR}/pufferfish.zip -d ${EXTERNAL_DIR}
mv ${EXTERNAL_DIR}/pufferfish-${SVER} ${EXTERNAL_DIR}/pufferfish

mkdir -p ${INSTALL_DIR}/include/pufferfish

cp ${EXTERNAL_DIR}/pufferfish/include/ProgOpts.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/BooPHF.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/SpinLock.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/Kmer.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/CanonicalKmer.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/string_view.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/CanonicalKmerIterator.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/PufferfishBaseIndex.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/PufferfishIndex.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/PufferfishSparseIndex.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/PufferfishLossyIndex.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/PufferfishTypes.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/rank9b.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/rank9sel.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/macros.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/select.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/Util.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/PairedAlignmentFormatter.hpp ${INSTALL_DIR}/include/pufferfish
#cp ${EXTERNAL_DIR}/pufferfish/include/SingleAlignmentFormatter.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/SelectiveAlignmentUtils.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/PuffAligner.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/MemCollector.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/MemChainer.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/CommonTypes.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/SAMWriter.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/PufferfishConfig.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/BulkChunk.hpp ${INSTALL_DIR}/include/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/include/BinWriter.hpp ${INSTALL_DIR}/include/pufferfish
cp -r ${EXTERNAL_DIR}/pufferfish/include/libdivide ${INSTALL_DIR}/include/pufferfish
cp -r ${EXTERNAL_DIR}/pufferfish/include/ksw2pp ${INSTALL_DIR}/include/pufferfish
cp -r ${EXTERNAL_DIR}/pufferfish/include/compact_vector ${INSTALL_DIR}/include/pufferfish
cp -r ${EXTERNAL_DIR}/pufferfish/include/metro ${INSTALL_DIR}/include/pufferfish
cp -r ${EXTERNAL_DIR}/pufferfish/include/chobo ${INSTALL_DIR}/include/pufferfish
cp -r ${EXTERNAL_DIR}/pufferfish/include/sparsepp ${INSTALL_DIR}/include/pufferfish
cp -r ${EXTERNAL_DIR}/pufferfish/include/simde ${INSTALL_DIR}/include/pufferfish
cp -r ${EXTERNAL_DIR}/pufferfish/include/tsl ${INSTALL_DIR}/include/pufferfish



mkdir -p ${INSTALL_DIR}/src/pufferfish
cp -r ${EXTERNAL_DIR}/pufferfish/src/metro ${INSTALL_DIR}/src/pufferfish
cp ${EXTERNAL_DIR}/pufferfish/src/rank9b.cpp ${INSTALL_DIR}/src/pufferfish


#mkdir -p ${INSTALL_DIR}/include/pufferfish
#mkdir -p ${INSTALL_DIR}/src/pufferfish
#
#rm ${EXTERNAL_DIR}/pufferfish/src/xxhash.c
#rm ${EXTERNAL_DIR}/pufferfish/include/xxhash.h
#
#cp -r ${EXTERNAL_DIR}/RapMap/src/*.c ${INSTALL_DIR}/src/rapmap
#cp -r ${EXTERNAL_DIR}/RapMap/src/*.cpp ${INSTALL_DIR}/src/rapmap
#cp -r ${EXTERNAL_DIR}/RapMap/src/metro ${INSTALL_DIR}/src/rapmap
#cp -r ${EXTERNAL_DIR}/RapMap/src/ksw2pp ${INSTALL_DIR}/src/rapmap
#cp -r ${EXTERNAL_DIR}/RapMap/include/tclap ${INSTALL_DIR}/include/rapmap
#cp -r ${EXTERNAL_DIR}/RapMap/include/*.h ${INSTALL_DIR}/include/rapmap
#cp -r ${EXTERNAL_DIR}/RapMap/include/*.hpp ${INSTALL_DIR}/include/rapmap
#cp -r ${EXTERNAL_DIR}/RapMap/include/sparsepp ${INSTALL_DIR}/include/rapmap
#cp -r ${EXTERNAL_DIR}/RapMap/include/digestpp ${INSTALL_DIR}/include/rapmap
#cp -r ${EXTERNAL_DIR}/RapMap/include/chobo ${INSTALL_DIR}/include/rapmap
#cp -r ${EXTERNAL_DIR}/RapMap/include/metro ${INSTALL_DIR}/include/rapmap
#cp -r ${EXTERNAL_DIR}/RapMap/include/ksw2pp ${INSTALL_DIR}/include/rapmap
#cp -r ${EXTERNAL_DIR}/RapMap/include/tsl ${INSTALL_DIR}/include/rapmap
#
###
## Remove some redundant files that might otherwise be duplicated
###
#rm ${INSTALL_DIR}/include/rapmap/FastxParser.hpp
#rm ${INSTALL_DIR}/include/rapmap/concurrentqueue.h
#rm ${INSTALL_DIR}/include/rapmap/FastxParserThreadUtils.hpp
#rm ${INSTALL_DIR}/src/rapmap/FastxParser.cpp

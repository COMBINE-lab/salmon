#!/bin/bash
source /hbb_exe/activate

set -e

CPATH=`pwd`
echo "[Drone test] current path : ${CPATH}"
echo "[Drone test] making quant test directory"

export PATH=/root/miniconda2/bin:$PATH
export PYTHONPATH=/root/miniconda2/lib/python2.7/site-packages

echo "[Drone test] run nextflow pipeline"

nextflow tests/test_quant.nf
# store the nextflow return value
nextflowret=$?
if [ $nextflowret -ne 0 ]; then
    echo "[Drone test]: nextflow pipeline test_quant.nf failed!"
    exit 1
fi

echo "[Drone test] echoing quants here"
grep "spearman" sim/*.json

mkdir -p "/mnt/ci_res/${DRONE_REPO}/${DRONE_COMMIT_BRANCH}/${DRONE_COMMIT_SHA}/sim"
cp sim/*.json "/mnt/ci_res/${DRONE_REPO}/${DRONE_COMMIT_BRANCH}/${DRONE_COMMIT_SHA}/sim/"

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

echo "[Drone test] echoing quants here"

grep "spearman" sim/*.json

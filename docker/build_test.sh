#! /bin/bash
SALMON_VERSION=1.10.0
TMPDIR=/mnt/scratch2/DELETE_ME_TEMP docker build --no-cache -t combinelab/salmon:${SALMON_VERSION} -t combinelab/salmon:latest .

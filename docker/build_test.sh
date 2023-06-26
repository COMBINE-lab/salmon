#! /bin/bash
SALMON_VERSION=1.10.2
TMPDIR=/mnt/scratch7/DELETE_ME_TEMP docker build --no-cache -t combinelab/salmon:${SALMON_VERSION} -t combinelab/salmon:latest .

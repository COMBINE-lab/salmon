#! /bin/bash
SALMON_VERSION=1.11.4
TMPDIR=/mnt/scratch7/DELETE_ME_TEMP docker build --no-cache -t combinelab/salmon:${SALMON_VERSION} -t combinelab/salmon:latest .

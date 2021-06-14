#! /bin/bash
SALMON_VERSION=1.5.1
docker build --no-cache -t combinelab/salmon:${SALMON_VERSION} -t combinelab/salmon:latest .

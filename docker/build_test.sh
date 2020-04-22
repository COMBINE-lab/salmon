#! /bin/bash
SALMON_VERSION=1.2.1
docker build --no-cache -t combinelab/salmon:${SALMON_VERSION} -t combinelab/salmon:latest .

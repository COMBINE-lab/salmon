#! /bin/bash
SALMON_VERSION=1.9.0
docker build --no-cache -t combinelab/salmon:${SALMON_VERSION} -t combinelab/salmon:latest .

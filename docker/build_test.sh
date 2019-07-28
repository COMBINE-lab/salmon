#! /bin/bash
SALMON_VERSION=0.15.0
docker build --no-cache -t combinelab/salmon:${SALMON_VERSION} -t combinelab/salmon:latest .

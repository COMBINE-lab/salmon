#! /bin/bash
SALMON_VERSION=0.12.1
docker build -t combinelab/salmon:${SALMON_VERSION} -t combinelab/salmon:latest .

#! /bin/bash
SALMON_VERSION=0.11.3
docker build -t combinelab/salmon:${SALMON_VERSION} -t combinelab/salmon:latest .

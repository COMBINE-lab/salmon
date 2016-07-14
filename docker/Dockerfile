# image: COMBINE-lab/salmon
# This dockerfile is based on the one created by
# Titus Brown (available at https://github.com/ctb/2015-docker-building/blob/master/salmon/Dockerfile)
FROM ubuntu:15.10
MAINTAINER salmon.maintainer@gmail.com

ENV PACKAGES git gcc make g++ cmake libboost-all-dev liblzma-dev libbz2-dev \
    ca-certificates zlib1g-dev curl unzip autoconf
ENV SALMON_VERSION v0.7.0-pre

# salmon binary will be installed in /home/salmon/bin/salmon

### don't modify things below here for version updates etc.

WORKDIR /home

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN curl -k -L https://github.com/COMBINE-lab/salmon/archive/v0.7.0-pre.tar.gz -o salmon-v0.7.0-pre.tar.gz && \
    tar xzf salmon-v0.7.0-pre.tar.gz && \
    cd salmon-0.7.0-pre && \
    mkdir build && \
    cd build && \
    cmake .. && make && make install

ENV PATH /home/salmon-0.7.0-pre/bin:${PATH}
FROM debian
MAINTAINER robp <rob.patro@cs.stonybrook.edu>
RUN apt-get update && apt-get -y upgrade \
    && apt-get -y install wget \
    && wget -O salmon.tar.gz https://github.com/COMBINE-lab/salmon/archive/v0.7.0-pre.tar.gz\
    && tar xvfz salmon.tar.gz \
    && apt-get clean \
    && rm -r /var/lib/apt/lists/*
WORKDIR /salmon-0.7.0-pre/bin
ENTRYPOINT ["./salmon"]
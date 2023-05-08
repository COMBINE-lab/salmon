# image: COMBINE-lab/pufferfish
# This dockerfile is based on the one created by
# Titus Brown (available at https://github.com/ctb/2015-docker-building/blob/master/salmon/Dockerfile)
#FROM insilicodb/bioconductor 
FROM ubuntu:16.04

ENV PACKAGES git gcc make g++ libboost-all-dev liblzma-dev libbz2-dev \
    ca-certificates zlib1g-dev curl unzip autoconf vim wget time bzip2
ENV PUFFERFISH_VERSION 0.9.0

# pufferfish executable will be available in /home/pufferfish/build/src
# salmon binary will be installed in /home/salmon/bin/salmon

### don't modify things below here for version updates etc.

WORKDIR /home

ADD https://cmake.org/files/v3.9/cmake-3.9.5-Linux-x86_64.sh /cmake-3.9.5-Linux-x86_64.sh
RUN mkdir /opt/cmake && \
	sh /cmake-3.9.5-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && \
	ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

RUN curl -o /usr/local/bin/jq http://stedolan.github.io/jq/download/linux64/jq && \
	chmod +x /usr/local/bin/jq

# python and required packages
RUN yes | apt-get install software-properties-common && \
	add-apt-repository ppa:deadsnakes/ppa && \
	apt-get update
RUN yes | apt-get install python3.6
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
	printf "\nyes\n\nyes\n" | bash Miniconda3-latest-Linux-x86_64.sh

RUN yes | /root/miniconda3/bin/conda install scipy && \
	yes | /root/miniconda3/bin/conda install pandas && \
	yes | /root/miniconda3/bin/conda install scikit-learn

# TwoPaCo
RUN wget https://github.com/01org/tbb/releases/download/2018_U4/tbb2018_20180411oss_lin.tgz && \
tar zxf tbb2018_20180411oss_lin.tgz && \
rm tbb2018_20180411oss_lin.tgz && \
cd tbb2018_20180411oss && \
cp -r include /usr/local && \
cp -r lib/intel64/gcc4.7/* /usr/local/lib && \
cd .. && \
rm -rf tbb2018_20180411oss

RUN git clone https://github.com/fataltes/TwoPaCo.git
RUN cd TwoPaCo && \
	git checkout pufferize && \
	mkdir build && \
	cd build && \
	cmake ../src && \
	make

RUN git clone https://github.com/COMBINE-lab/pufferfish.git

# sdsl-lite
RUN	cd pufferfish && \
	git clone https://github.com/simongog/sdsl-lite.git && \
	cd sdsl-lite && \ 
	./install.sh ../
	
# Pufferfish
RUN cd pufferfish && \	
	git checkout develop && \
	#git checkout b7f383519d6fcad1d6392d7502e1bef23c0c36ed && \
	mkdir build && \
	cd build && \
	cmake ../ && \
	make

# Salmon
#RUN git clone https://github.com/fataltes/salmon.git && \
#	cd salmon && \
#	git checkout salpuf && \
#	mkdir build && \
#	cd build && \
#	cmake -DFETCH_BOOST=TRUE ../ && \
#	make && \
#	make install

RUN ldconfig

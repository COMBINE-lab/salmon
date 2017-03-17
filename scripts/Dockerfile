FROM combinelab/holy-build-box-64:latest_working

RUN set -x

# Install things we need
RUN yum install -y --quiet wget
RUN wget http://download.fedoraproject.org/pub/epel/5/x86_64/epel-release-5-4.noarch.rpm
RUN rpm -i --quiet epel-release-5-4.noarch.rpm
#yum install -y --quiet git
RUN yum install -y --quiet unzip
RUN yum install -y --quiet bzip2-devel.x86_64
RUN yum install -y --quiet xz-devel.x86_64
RUN yum install -y --quiet git

RUN wget http://downloads.sourceforge.net/project/boost/boost/1.59.0/boost_1_59_0.tar.gz 
RUN tar xzf boost_1_59_0.tar.gz
WORKDIR "/boost_1_59_0"
RUN source /hbb_exe/activate && ./bootstrap.sh --prefix=/usr --with-libraries=iostreams,atomic,chrono,container,date_time,exception,filesystem,graph,graph_parallel,math,program_options,system,thread,timer,serialization
RUN source /hbb_exe/activate && ./b2 -d0 -j10 cxxflags=-std=c++11 link=static install
WORKDIR "/"
RUN rm boost_1_59_0.tar.gz
RUN rm -fr "/boost_1_59_0"
RUN wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
RUN bash Miniconda2-latest-Linux-x86_64.sh -b
RUN source /hbb_exe/activate && PYTHONPATH=/root/miniconda2/lib/python2.7/site-packages PATH=/root/miniconda2/bin:$PATH pip install pandas scipy numpy 
#matplotlib seaborn

# java
RUN wget --no-cookies --no-check-certificate --header "Cookie: gpw_e24=http%3A%2F%2Fwww.oracle.com%2F; oraclelicense=accept-securebackup-cookie" "http://download.oracle.com/otn-pub/java/jdk/8u60-b27/jre-8u60-linux-x64.rpm" -O jre-8u60-linux-x64.rpm
RUN yum localinstall --nogpgcheck -y --quiet jre-8u60-linux-x64.rpm
RUN rm jre-8u60-linux-x64.rpm

# and nextflow
RUN curl -fsSL get.nextflow.io | bash
RUN mv nextflow /usr/local/bin/

RUN yum install -y --quiet shasum 

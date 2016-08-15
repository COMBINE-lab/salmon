#!/bin/bash
set -e

branch=$1
version=$2

echo "Building salmon [branch = ${branch}]. Tagging version as ${version}"

# Activate Holy Build Box environment.
source /hbb_exe/activate

set -x

# Install things we need
yum install -y --quiet wget
wget http://download.fedoraproject.org/pub/epel/5/x86_64/epel-release-5-4.noarch.rpm
rpm -i --quiet epel-release-5-4.noarch.rpm
#yum install -y --quiet git
yum install -y --quiet unzip
yum install -y --quiet bzip2-devel.x86_64
yum install -y --quiet xz-devel.x86_64

curl -k -L https://github.com/COMBINE-lab/salmon/archive/${branch}.zip -o ${branch}.zip
unzip ${branch}.zip
mv salmon-${branch} salmon
cd salmon
mkdir build
cd build
cmake -DFETCH_BOOST=TRUE ..
make
make install
make test
cd ../scripts
bash make-release.sh -v ${version} -n linux_x86_64
cd ../RELEASES
cp *.tar.gz /io/

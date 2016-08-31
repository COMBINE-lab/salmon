#~/bin/bash

echo "building pre-compiled linux release for Salmon $1"

docker run -t -i --rm   -v `pwd`:/io   phusion/holy-build-box-64:latest   bash /io/compile.sh develop $1

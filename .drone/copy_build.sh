#!/bin/bash
source /hbb_exe/activate

set -e

cd ../scripts
bash make-release.sh -v latest -n linux_x86_64
cd ../RELEASES
mkdir -p "/mnt/data/${DRONE_REPO}/${DRONE_REPO_BRANCH}/build"
cp *.tar.gz "/mnt/data/${DRONE_REPO}/${DRONE_REPO_BRANCH}/build/"

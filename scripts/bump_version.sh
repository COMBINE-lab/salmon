#!/usr/bin/env bash

HERE="$( cd "$(dirname "$0")" ; pwd -P )"

major_v=`${HERE}/../build/src/salmon -v | cut -d ' ' -f 2 | cut -d '.' -f 1`
minor_v=`${HERE}/../build/src/salmon -v | cut -d ' ' -f 2 | cut -d '.' -f 2`
patch_v=`${HERE}/../build/src/salmon -v | cut -d ' ' -f 2 | cut -d '.' -f 3`

echo "VERSION : ${major_v}.${minor_v}.${patch_v}"

# update docker file
awk -v majv=${major_v} -v minv=${minor_v} -v patchv=${patch_v} '{ if ($0 ~ /ENV SALMON_VERSION/) { print "ENV SALMON_VERSION " majv"."minv"."patchv; } else { print $0; }}' \
    ${HERE}/../docker/Dockerfile > ${HERE}/../docker/Dockerfile.new && mv ${HERE}/../docker/Dockerfile.new ${HERE}/../docker/Dockerfile

# update docker build script
awk -v majv=${major_v} -v minv=${minor_v} -v patchv=${patch_v} '{ if ($0 ~ /SALMON_VERSION=/) { print "SALMON_VERSION="majv"."minv"."patchv; } else { print $0; }}' \
    ${HERE}/../docker/build_test.sh > ${HERE}/../docker/build_test.sh.new && mv ${HERE}/../docker/build_test.sh.new ${HERE}/../docker/build_test.sh

# update version file (which feeds cmake)
echo -e "VERSION_MAJOR ${major_v}\nVERSION_MINOR ${minor_v}\nVERSION_PATCH ${patch_v}" > ${HERE}/../current_version.txt

# update conf.py
awk -v majv=${major_v} -v minv=${minor_v} -v patchv=${patch_v} \
    '{ if ($0 ~ /version = /) { print "version = '\''" majv"."minv"'\''"; } else if ($0 ~ /release = /) { print "release = '\''" majv"."minv"."patchv"'\''"; } else { print $0; }}' \
    ${HERE}/../doc/source/conf.py > ${HERE}/../doc/source/conf.py.new && mv ${HERE}/../doc/source/conf.py.new ${HERE}/../doc/source/conf.py

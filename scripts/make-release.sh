#!/bin/bash

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

host=
version=

while getopts "v:n:" opt; do
  case $opt in
    n)
      echo "Host is $OPTARG" >&2
      host=$OPTARG
      ;;
    v)
      echo "Version is $OPTARG" >&2
      version=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

echo -e "Preparing binary release\n=====================\n"
echo -e "Version = ${version}"
echo -e "Host = ${host}"

# create the binary directory 
betaname=salmon-${version}_${host}
mkdir ${DIR}/../RELEASES
mkdir ${DIR}/../RELEASES/${betaname}
mkdir ${DIR}/../RELEASES/${betaname}/bin
mkdir ${DIR}/../RELEASES/${betaname}/lib

# copy over the executable and Intel TBB libraries
echo -e "Copying over the binary and Intel TBB libraries\n"
cp ${DIR}/../bin/salmon ${DIR}/../RELEASES/${betaname}/bin/
cp ${DIR}/../lib/libtbb* ${DIR}/../RELEASES/${betaname}/lib/

# copy over the test data
cp ${DIR}/../sample_data.tgz ${DIR}/../RELEASES/${betaname}/

# copy other dependencies (shared libraries)
echo -e "Copying over other shared library dependencies\n"
bash ${DIR}/../scripts/cpld.bash ${DIR}/../bin/salmon ${DIR}/../RELEASES/${betaname}/lib/
echo -e "Removing dangerous dependencies\n"
rm ${DIR}/../RELEASES/${betaname}/lib/librt.so.1
rm ${DIR}/../RELEASES/${betaname}/lib/libc.so.6
rm ${DIR}/../RELEASES/${betaname}/lib/ld-linux-x86-64.so.2
rm ${DIR}/../RELEASES/${betaname}/lib/libdl.so.2
rm ${DIR}/../RELEASES/${betaname}/lib/libpthread*.so.*

# now make the tarball
echo -e "Making the tarball\n"
cd ${DIR}/../RELEASES
tar czvf ${betaname}.tar.gz ${betaname}

echo -e "Done making release!"

#echo -e "Pushing the tarball to GitHub\n"
## Since it's currently unclear to me how to overwrite an asset via the GitHub
## API, the following code deletes the old asset, and uploads the new one in its place
#
## Get the previous asset id of the tarball
#echo -e "Getting previous asset ID\n"
#ASSETID=`curl -s -X GET https://api.github.com/repos/COMBINE-lab/salmon/releases/1263754/assets | grep "\"id" | head -1 | awk '{gsub(/,$/,""); print $2}'`
#
## Delete the previous tarball
#echo -e "Deleting previous asset\n"
#curl -X DELETE -H "Authorization: token ${SALMON_PUSH_KEY}" https://api.github.com/repos/COMBINE-lab/salmon/releases/assets/$ASSETID
#
## Upload the new tarball
#echo -e "Uploading new asset\n"
#curl -X POST --data-binary "@SalmonBeta-latest_ubuntu-12.04.tar.gz" https://uploads.github.com/repos/COMBINE-lab/salmon/releases/1263754/assets?name=SalmonBeta-latest_ubuntu-12.04.tar.gz --header "Content-Type:application/gzip" -H "Authorization: token ${SALMON_PUSH_KEY}"
#echo -e "Done!\n"

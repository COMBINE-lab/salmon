
major_v=`../build/src/salmon -v | cut -d ' ' -f 2 | cut -d '.' -f 1`
minor_v=`../build/src/salmon -v | cut -d ' ' -f 2 | cut -d '.' -f 2`
patch_v=`../build/src/salmon -v | cut -d ' ' -f 2 | cut -d '.' -f 3`

echo "VERSION : ${major_v}.${minor_v}.${patch_v}"

awk -v majv=${major_v} -v minv=${minor_v} -v patchv=${patch_v} '{ if ($0 ~ /ENV SALMON_VERSION/) { print "ENV SALMON_VERSION " majv"."minv"."patchv; } else { print $0; }}' ../docker/Dockerfile > ../docker/Dockerfile.new

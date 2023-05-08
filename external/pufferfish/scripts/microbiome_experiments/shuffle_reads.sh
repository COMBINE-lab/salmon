r1=$1
r2=$2

cmd="paste -d \$'\10' <(awk '{printf(\"%s%s\",\$0,(NR%4==0)?\"\n\":\"\0\")}' ${r1}) <(awk '{printf(\"%s%s\",\$0,(NR%4==0)?\"\n\":\"\0\")}' ${r2}) | shuf | awk -v FS=\$'\\10' '{ print \$1 > \"${r1}.shuf\" ; print \$2 > \"${r2}.shuf\" }'"
echo ${cmd}
eval ${cmd}

cmd="cat ${r1}.shuf | tr \"\0\" \"\n\" > tmp.fastq"
echo ${cmd}
eval ${cmd}
mv tmp.fastq ${r1}
rm ${r1}.shuf

cmd="cat ${r2}.shuf | tr \"\0\" \"\n\" > tmp.fastq"
echo ${cmd}
eval ${cmd}
mv tmp.fastq ${r2}
rm ${r2}.shuf

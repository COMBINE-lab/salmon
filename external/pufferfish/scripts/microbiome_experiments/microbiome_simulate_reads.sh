
### User only needs to set these variables
root="/mnt/scratch3/SubmitPuffalign"
minReadCnt=20
finalCnt=1000000
###

if [[ $# -eq 0 ]]; then
    echo "Need to have the address to where the sample real paired-end read files  are relative to the root address"
	exit
fi

dataDir="${root}/$1"
r1="${dataDir}/S1R1.fq"
r2="${dataDir}/S1R2.fq"

if [[ ! -d ${dataDir} ]]; then
	echo "${dataDir} does not exist"
	exit
fi

if [[ ! -f "${r1}" || ! -f "${r2}" ]]; then
	echo "${r1} or ${r2} or both do not exist"
	echo "make sure you've renamed your read files to S1R1.fq and S1R2.fq or change the read names in this script accordingly."
	exit
fi

simulatedReads="${dataDir}/simulatedReads"
realReads="${dataDir}/realReads_bowtie2Aligned"
mkdir -p ${dataDir}/truth

### STEP1
echo "STEP1: Run Bowtie2"
cmd="/usr/bin/time ${root}/tools/bowtie2-2.3.5.1-linux-x86_64/bowtie2 -k 1 --maxins 1000 --gbar 1 --no-discordant --no-mixed  -x ${root}/indices/bowtie2indices/microbiome.bowtie2index/bw2 -1 ${r1} -2 ${r2} -S ${realReads}.sam -p 16 --ignore-quals > ${realReads}.log 2>&1"
echo ${cmd}
eval ${cmd}

### STEP2
echo ""
echo "STEP2: Scale count up to 1M reads"
currCnt=$(samtools view -F 4 ${realReads}.sam | cut -f 3 | sort | uniq -c | awk -v minReadCnt="${minReadCnt}" 'BEGIN {s=0} {if ($1 >= minReadCnt*2) {s+=$1}} END {print s}')
scale=1
if (( ${currCnt} < ${finalCnt} )); then
	scale=$(( (finalCnt/currCnt) + 1 ))
fi
echo "current Cnt=${currCnt} and scale=${scale}"


### STEP3
echo ""
echo "STEP3: Generate Profile"
cmd="samtools view -F 4 ${realReads}.sam | cut -f 3 | sort | uniq -c | awk '{if (\$1 >= ${minReadCnt}*2) {print \$1*${scale}\" \"\$2}}' | sort -nr > ${simulatedReads}.profile"
echo ${cmd}
eval ${cmd}


### STEP4
echo ""
echo "STEP4: Check if all the reference files exist"
foundAll=1
totalRefs=0
while read l;do
	count=$(echo ${l} | cut -d' ' -f1)
	ref=$(echo ${l} | cut -d' ' -f2 | sed 's/|/_/g; s/kraken:taxid_//g')
	if [ ! -f "${root}/microbiome/bacteria_virus/library/bacteria/splitted_selected/splitted/$ref.fna" ]; then
		echo "${root}/microbiome/bacteria_virus/library/bacteria/splitted_selected/splitted/$ref.fna does not exist"
		foundAll=0
	fi
	totalRefs=$(( totalRefs + 1 ))
done < ${simulatedReads}.profile

if (( $foundAll == 0 )); then
	exit
fi
echo "Yup! They do!"
echo "We have all the ${totalRefs} references present"


### STEP5
refCntr=0
echo ""
echo ""
echo "STEP5: Create mason simulated reads for each of the references and combine all in one file"
echo ""
while read l;do
	count=$(echo ${l} | cut -d' ' -f1)
	ref=$(echo ${l} | cut -d' ' -f2 | sed 's/|/_/g; s/kraken:taxid_//g')
	refCntr=$(( refCntr + 1 ))
	echo "${refCntr} ==> count: ${count} ref: ${ref}"
	cmd="${root}/tools/mason2-2.0.9-Linux-x86_64/bin/mason_simulator --embed-read-info -oa ${dataDir}/truth/${ref}.bam -ir ${root}/microbiome/bacteria_virus/library/bacteria/splitted_selected/splitted/$ref.fna -n $count -o ${dataDir}/tmp_1.fastq -or ${dataDir}/tmp_2.fastq --read-name-prefix $ref"
#echo ${cmd}
	eval ${cmd}
	cmd="cat ${dataDir}/tmp_1.fastq >> ${simulatedReads}_1.fastq"
#echo ${cmd}
	eval ${cmd}
	cmd="cat ${dataDir}/tmp_2.fastq >> ${simulatedReads}_2.fastq"
#echo ${cmd}
	eval ${cmd}
	rm ${dataDir}/tmp_1.fastq
	rm ${dataDir}/tmp_2.fastq
done < ${simulatedReads}.profile


### STEP6
echo ""
echo ""
echo "STEP6: Shuffle simulated reads"
bash shuffle_reads.sh ${simulatedReads}_1.fastq ${simulatedReads}_2.fastq
echo ""
echo "DONE SIMULATING READS"

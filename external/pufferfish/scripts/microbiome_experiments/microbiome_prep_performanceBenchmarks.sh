#!/bin/bash
if [[ ! $# -eq 1 ]]; then
	echo "ERROR! Need the type of the experiment. Accepted inputs: \"bulk\" and \"singleStrain\""
	exit
fi

scriptDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

source ${scriptDir}/shared_setting.sh

if [[ $1 == "bulk" ]]; then
	echo "---------------------> Running for BULK"
	source ${scriptDir}/bulk_setting.sh
elif [[ $1 == "singleStrain" ]]; then
	echo "---------------------> Running for SINGLE_STRAIN"
	source ${scriptDir}/singleStrain_setting.sh
else
	echo "ERROR! Wrong argument passed. The only valid arguments are \"bulk\" and \"singleStrain\""
	exit
fi

options=( 'primary' 'maxNumHits20' 'maxNumHits200' )
if [[ ${experimentType} == "bulk" ]]; then
	outfile="${root}/results/bulk_performance.res"
elif [[ ${experimentType} == "singleStrain" ]]; then
	outfile="${root}/results/singleStrain_performance.res"
else
	echo "ERROR! Unknown experiment type: ${experimentType}. Valid values are bulk or singleStrain"
	exit
fi

if [ -f ${outfile} ]; then
	rm ${outfile}
fi

for option in ${options[@]}; do
	for dataType in ${samples[@]}; do
		for wc in 'warm' 'cold'; do
			echo "${dataType} ${option} ${wc}"
			v=$(grep 'maxresident' ${root}/alignments/puffaligner/${dataType}.pe.concordant.${option}_${wc}.log)
			totalTime=$(echo ${v} | tr ' ' '\n' | grep 'elapsed' | sed 's/elapsed//')
			mem=$(echo ${v} | tr ' ' '\n' | grep 'maxresident' | sed 's/maxresident)k//')
			alignTime=$(grep 'Elapsed time' ${root}/alignments/puffaligner/${dataType}.pe.concordant.${option}_${wc}.log | sed 's/Elapsed time: //; s/s//;')
			echo "${wc} ${option} ${dataType} puffaligner ${mem} ${totalTime} ${alignTime}" >> ${outfile}

			v=$(grep 'maxresident' ${root}/alignments/bowtie2/${dataType}.pe.concordant.${option}_${wc}.log)
			totalTime=$(echo ${v} | tr ' ' '\n' | grep 'elapsed' | sed 's/elapsed//')
			mem=$(echo ${v} | tr ' ' '\n' | grep 'maxresident' | sed 's/maxresident)k//')
			alignTime=$(grep 'Multiseed full-index search' ${root}/alignments/bowtie2/${dataType}.pe.concordant.${option}_${wc}.log | sed 's/Multiseed full-index search: //')
			echo "${wc} ${option} ${dataType} bowtie2 ${mem} ${totalTime} ${alignTime}" >> ${outfile}

			if [[ ${runStar} == 1 ]]; then
				v=$(grep 'maxresident' ${root}/alignments/star/${dataType}.pe.concordant.${option}/star_${wc}.log)
				totalTime=$(echo ${v} | tr ' ' '\n' | grep 'elapsed' | sed 's/elapsed//')
				mem=$(echo ${v} | tr ' ' '\n' | grep 'maxresident' | sed 's/maxresident)k//')
				alignTime1=$(grep 'started mapping' ${root}/alignments/star/${dataType}.pe.concordant.${option}/star_${wc}.log | cut -d' '  -f3)
				alignTime2=$(grep 'finished mapping' ${root}/alignments/star/${dataType}.pe.concordant.${option}/star_${wc}.log | cut -d' '  -f3)
				alignTime=${alignTime1}-${alignTime2}
				echo "${wc} ${option} ${dataType} star ${mem} ${totalTime} ${alignTime}" >> ${outfile}
			fi
			
			v=$(grep 'maxresident' ${root}/alignments/deBGA/${dataType}.pe.concordant.${option}_${wc}.log)
			totalTime=$(echo ${v} | tr ' ' '\n' | grep 'elapsed' | sed 's/elapsed//')
			mem=$(echo ${v} | tr ' ' '\n' | grep 'maxresident' | sed 's/maxresident)k//')
			alignTime=$(grep 'seconds is used in mapping' ${root}/alignments/deBGA/${dataType}.pe.concordant.${option}_${wc}.log | sed 's/ seconds is used in mapping//')
			echo "${wc} ${option} ${dataType} debga ${mem} ${totalTime} ${alignTime}" >> ${outfile}

		done
	done
done

#!/bin/bash

#### NOTE BEFORE RUNNING THIS SCRIPT ####
###	1. make sure you have built all the corresponding indices for bulk and singleStrain experiments
### 2. make sure you have correctly set the values for variables in shared_setting.sh , bulk_setting.sh , and singleStrain_setting.sh
### 3. make sure all the used scripts are in the same directory as the current bash script
### 4. ENJOY!

scriptRoot="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "${scriptRoot}"

source ${scriptRoot}/shared_setting.sh


echo ""
echo ""
echo "-----------------------------> SINGLE_STRAIN experiment <-----------------------------"
source ${scriptRoot}/singleStrain_setting.sh
echo "Align and Quantify Samples"
cmd="bash ${scriptRoot}/microbiome_align_quantify.sh singleStrain"; echo ${cmd}; eval ${cmd}; echo "";

echo "Write down accuracy results"
cmd="python ${scriptRoot}/clean_results.py quant ${root}/quants/singleStrain.res ${root}/results/singleStrain_accuracy.xlsx"; echo ${cmd}; eval ${cmd}; echo "";

echo "Write down performance results"
cmd="bash ${scriptRoot}/microbiome_prep_performanceBenchmarks.sh singleStrain"; echo ${cmd}; eval ${cmd}; echo "";
cmd="python ${scriptRoot}/clean_results.py perf ${root}/results/singleStrain_performance.res ${root}/results/singleStrain_performance.xlsx"; echo ${cmd}; eval ${cmd}; echo "";

echo ""
echo ""
echo "-----------------------------> BULK experiment <-----------------------------"
source ${scriptRoot}/bulk_setting.sh
echo "Align and Quantify Samples"
cmd="bash ${scriptRoot}/microbiome_align_quantify.sh bulk"; echo ${cmd}; eval ${cmd}; echo "";

echo "Evaluate and write down the accuracy results"
sampleList=""
for sample in ${samples[@]}; do
	sampleList+="${sample},"
done
cmd="python ${scriptRoot}/evaluate_bulk_results.py ${root} ${readDir_relative2root} ${sampleList} ${readName}"; echo ${cmd}; eval ${cmd}; echo "";

echo "Write down performance results"
cmd="bash ${scriptRoot}/microbiome_prep_performanceBenchmarks.sh bulk"; echo ${cmd}; eval ${cmd}; echo "";
cmd="python ${scriptRoot}/clean_results.py perf ${root}/results/bulk_performance.res ${root}/results/bulk_performance.xlsx"; echo ${cmd}; eval ${cmd}; echo "";

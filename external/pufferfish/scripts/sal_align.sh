
config='config.json'
if [[ $# -eq 1 ]]; then
	config=$1
fi

pufferfish=$(jq -r '.pufferfish' ${config})/build/src/pufferfish
salmon=$(jq -r '.salmon' ${config})/build/src/salmon
salidx=$(jq -r '.salidx' ${config})
fasta=$(jq -r '.fasta' ${config})
ratio1=$(jq -r '.ratio' ${config})
read1=$(jq -r '.read1' ${config})
read2=$(jq -r '.read2' ${config})
priorr=$(jq -r '.prior' ${config})
outdir=$(jq -r '.outDir' ${config})
outfile=$(jq -r '.outFile' ${config})
#outfile=${outfile}_sr${ratio1}
truthfile=$(jq -r '.truthFile' ${config})
taxatruthfile=$(jq -r '.taxaTruthFile' ${config})


echo "salmon normal EM .."
echo "/usr/bin/time ${salmon} quant -p 16 -la -i ${salidx} -1 ${read1} -2 ${read2}  -o ${outdir}/results/salmon/${outfile} --meta"
/usr/bin/time ${salmon} quant -p 16 -la -i ${salidx} -1 ${read1} -2 ${read2}  -o ${outdir}/results/salmon/${outfile} --meta

echo "reference level metric results:"
echo "python ${outdir}/scripts/validate_results.py -truer ${truthfile} -est ${outdir}/results/salmon/${outfile}/quant.sf -truthType microbiome"
~/miniconda3/bin/python ${outdir}/scripts/validate_results.py -truer ${truthfile} -est ${outdir}/results/salmon/${outfile}/quant.sf -truthType microbiome

echo "taxonomy metric results:"
mkdir -p ${outdir}/results/script_results/${outfile}_sal                             
echo "python ${outdir}/scripts/taxonomicRanks.py -qf ${outdir}/results/salmon/${outfile}/quant.sf -lf ${outdir}/files/rankedlineage.dmp -od ${outdir}/results/script_results/${outfile}_sal -dt salmon"
~/miniconda3/bin/python ${outdir}/scripts/taxonomicRanks.py -qf ${outdir}/results/salmon/${outfile}/quant.sf -lf ${outdir}/files/rankedlineage.dmp -od ${outdir}/results/script_results/${outfile}_sal -dt salmon

echo "python ${outdir}/scripts/validate_results.py -truer ${taxatruthfile} -est ${outdir}/results/script_results/${outfile}_sal/out.tab -truthType taxonomy"
~/miniconda3/bin/python ${outdir}/scripts/validate_results.py -truer ${taxatruthfile} -est ${outdir}/results/script_results/${outfile}_sal/out.tab -truthType taxonomy


#echo "salmon vb .."

#/usr/bin/time ${salmon} quant -p 16 -la -i ${salidx} -1 ${read1} -2 ${read2} -o ${outdir}/results/salmon/${outfile}_vb${priorr} --meta --useVBOpt --vbPrior 1e-${priorr}

#mkdir -p ${outdir}/results/script_results/${outfile}_vb${priorr}_sal
#python ${outdir}/scripts/taxonomicRanks.py -qf ${outdir}/results/salmon/${outfile}_vb${priorr}/quant.sf -lf ${outdir}/files/rankedlineage.dmp -od ${outdir}/results/script_results/${outfile}_vb${priorr}_sal
#python ${outdir}/scripts/calc_metrics.py -tr ${outdir}/files/true_counts.tab -s ${outdir}/results/script_results/${outfile}_vb${priorr}_sal/out.tab


#!/usr/bin/env nextflow

params.salmon = '/drone/src/github.com/COMBINE-lab/salmon/bin/salmon'
ref = file('/mnt/data/simulated/sim/Homo_sapiens.GRCh37.75.cdna.pc.fa')
truthpath = file('/mnt/data/simulated/sim/truth')
basepath = Channel.from('/mnt/data/simulated/sim/out/out')
scriptdir = file('/drone/src/github.com/COMBINE-lab/salmon/scripts')
resdir = './'
conds = ['A', 'B']
samples = [1, 2]

process buildIndex {
        cpus 2
        input:
        file ref

        output:
        file nfindex into index

        """
        ${params.salmon} index -t $ref -i nfindex
        """
}

process quantSim {
  cpus 16

  input:
  file index
  val basepath 
  each cond from conds
  each sample from samples

  output:
  file "sim_quants/${cond}_${sample}" into simqs

	script:
	if (cond == 'A')
	    """
      ${params.salmon} quant -p 16 -i ${index} -l A -1 ${basepath}_${sample}/sample_01_1_shuffled.fa.gz -2 ${basepath}_${sample}/sample_01_2_shuffled.fa.gz -o sim_quants/${cond}_${sample} --useVBOpt --validateMappings --minScoreFraction 0.95  --rangeFactorizationBins 4 --vbPrior 1e-2 --perTranscriptPrior
  	    """
	else
	    """
      ${params.salmon} quant -p 16 -i ${index} -l A -1 ${basepath}_${sample}/sample_02_1_shuffled.fa.gz -2 ${basepath}_${sample}/sample_02_2_shuffled.fa.gz -o sim_quants/${cond}_${sample} --useVBOpt --validateMappings --minScoreFraction 0.95 --rangeFactorizationBins 4 --vbPrior 1e-2 --perTranscriptPrior
	    """
}

process evalQuants {

  publishDir "$resdir"

	input:
  val truthpath
	file flist from simqs.toList()
  	
  output:
  file sim into simres

	script:
	"""
	for f in ${flist}
	do
    if [[ \$f == A* ]]
    then
      python ${scriptdir}/test_sim_corr.py --sim "\${f}/quant.sf" --est "${truthpath}/truthA.tpm" --out "sim/\${f}_res.json"
    else
      python ${scriptdir}/test_sim_corr.py --sim "\${f}/quant.sf" --est "${truthpath}/truthB.tpm" --out "sim/\${f}_res.json"
    fi
	done
	"""
}

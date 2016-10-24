#!/bin/bash
for fn in data/DRR0161{25..40};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
sailfish quant -i athal_index -l A \
         -1 ${fn}/${samp}_1.fastq.gz \
         -2 ${samp}_2.fastq.gz \
		 -p 8 -o quants/${samp}_quant
done 

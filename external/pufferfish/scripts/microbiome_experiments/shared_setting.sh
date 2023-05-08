#!/bin/bash

root="/mnt/scratch3/SubmitPuffalign"

##### Alignment Tools #####
puff="/home/mohsen/Pufferfish-repo/pufferfish/build/src/pufferfish"
bt2="${root}/tools/bowtie2-2.3.5.1-linux-x86_64/bowtie2"
star="${root}/tools/STAR-2.7.3a/bin/Linux_x86_64/STAR"	
debga="${root}/tools/deBGA/deBGA"
tools=( ${puff} ${bt2} ${star} ${debga} )
## Salmon as quantifier ##
salmon="${root}/tools/salmon-latest_linux_x86_64/bin/salmon"

##### variable Options #####
options=( 'primary' 'maxNumHits20' 'maxNumHits200' )
puffOpts=( '--bestStrata --primaryAlignment' '--maxNumHits 20' '--maxNumHits 200' )
bw2Opts=( '-k 1' '-k 20' '-k 200' )
starOpts=( '--outSAMmultNmax 1' '--outSAMmultNmax 20' '--outSAMmultNmax 200' )
debgaOpts=( '-o 1 -x 1' '-o 20 -x 20' '-o 200 -x 200' )

##### Output directories #####
mkdir -p ${root}/alignments/puffaligner
mkdir -p ${root}/alignments/bowtie2
mkdir -p ${root}/alignments/star
mkdir -p ${root}/alignments/deBGA

mkdir -p ${root}/quants/salmon
mkdir -p ${root}/results

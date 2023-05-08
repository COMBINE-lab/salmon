#!/bin/bash

##### root directory where all the data exists, all the experiments will be run and all results will be stored in some subdir of it. #####
root="/mnt/scratch3/SubmitPuffalign"

##### bulk vs sameSpecies #####
experimentType="singleStrain"

samples=( 'covid19' 'sars' 'bat2008' )

##### reference file that all the indices have been constructed over #####
referenceFile_relative2root="microbiome/bacteria_virus/virus2/library/viral/library.fna"

##### viral indices #####
puffIdx="${root}/indices/pufferfishindices/viralCG.puffindex"
bt2Idx="${root}/indices/bowtie2indices/viralCG.bowtie2index/bw2"
starIdx="${root}/indices/starindices/viralCG.starindex.fat6/"
debgaIdx="${root}/indices/deBGAindices/viralCG.deBGAindex/"

#### sample read addresses ####
readDir_relative2root="data/microbiome_simulation/singleStrainSample/coronas"
readName="simulatedReads"
#### The two read files are assumed to be the following:
########	r1="${root}/${readDir_relative2root}/${sample}/${readName}_1.fastq"
########	r2="${root}/${readDir_relative2root}/${sample}/${readName}_2.fastq"



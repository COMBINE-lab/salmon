#!/bin/bash

##### root directory where all the data exists, all the experiments will be run and all results will be stored in some subdir of it. #####
root="/mnt/scratch3/SubmitPuffalign"

##### bulk vs sameSpecies #####
experimentType="bulk"

samples=( 'jiandongPeninsula_SRR11283975' 'uzonCaldera_SRR11496426' 'soil_SRR1094822' )

##### reference file that all the indices have been constructed over #####
referenceFile_relative2root="microbiome/bacteria_virus/library/bacteria/selected_complete_genome.fna"

##### bulk indices #####
puffIdx="${root}/indices/pufferfishindices/microbiome_k31.resubmission.puffindex.sparse/"
bt2Idx="${root}/indices/bowtie2indices/microbiome.bowtie2index/bw2"
starIdx="${root}/indices/starindices/microbiome.starindex.fat9/"
debgaIdx="${root}/indices/deBGAindices/microbiome.idx/"

#### sample read addresses ####
readDir_relative2root="data/microbiome_simulation/bulkSample"
readName="simulatedReads"
#### The two read files are assumed to be the following:
########	r1="${root}/${readDir_relative2root}/${sample}/${readName}_1.fastq"
########	r2="${root}/${readDir_relative2root}/${sample}/${readName}_2.fastq"
#### The truth file for the corresponding read (which is already generated via the script for simulating the reads) is assumed to be called:
#######		${root}/${readDir_relative2root}/${sample}/${readName}.profile


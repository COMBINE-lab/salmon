#!/bin/bash
mkdir data
cd data
for i in `seq 25 40`; 
do 
  mkdir DRR0161${i}; 
  cd DRR0161${i}; 
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR0161${i}/DRR0161${i}_1.fastq.gz; 
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016/DRR0161${i}/DRR0161${i}_2.fastq.gz; 
  cd ..; 
done
cd .. 

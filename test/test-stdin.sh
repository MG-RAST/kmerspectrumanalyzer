#!/bin/bash
set -x #echo on

rm *.21 

# download the datasets if needed
make NC_017625.fna SRR001665_1.fastq SRR001665_2.fastq.gz

# build the datasets 
cat NC_017625.fna NC_017625.fna > DOUBLE.fna

countkmer21.sh DOUBLE.fna          ; mv DOUBLE.fna.21        DOUBLE-single.21
countkmer21.sh NC_017625.fna       ; mv NC_017625.fna.21     SINGLE-single.21   
countkmer21.sh SRR001665_1.fastq   ; mv SRR001665_1.fastq.21 SRR_1-single.21
countkmer21.sh SRR001665_2.fastq.gz; mv SRR001665_2.fastq.21 SRR_2-single.21

#cat NC_017625.fna                | countkmer21.sh > SINGLE-stdin.21   # FAILS fasta not supported
#cat NC_017625.fna NC_017625.fna  | countkmer21.sh > DOUBLE-stdin.21   # FAILS fasta not supported
cat SRR001665_1.fastq            | countkmer21.sh > SRR_1-stdin.21 

countkmer21.sh NC_017625.fna DOUBLE.fna SRR001665_1.fastq SRR001665_2.fastq.gz  

mv NC_017625.fna.21      SINGLE-multiple.21
mv DOUBLE.fna.21         DOUBLE-multiple.21
mv SRR001665_1.fastq.21  SRR_1-multiple.21 
mv SRR001665_2.fastq.21  SRR_2-multiple.21 



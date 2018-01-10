#!/bin/bash
set -x #echo on

# This script tests the command-line counting functionality for several sequence formats.
# mostly, this checks that countkmer21.sh and kmer-tool2 give output when expected

# download the datasets if needed
make NC_017625.fna SRR033548.fastq  SRR033548.fastq.gz 

ln  SRR033548.fastq.gz  SRR033548A.fastq.gz   # need different names, but don't mind having the same data

# build the datasets 
cat NC_017625.fna NC_017625.fna  NC_017625.fna NC_017625.fna  NC_017625.fna NC_017625.fna NC_017625.fna NC_017625.fna NC_017625.fna NC_017625.fna > GENOME10X.fna
# Generate a 10x duplicate of this dataset 
for i in $(seq 10); do cat SRR001665_1.fastq  >> BIG.fastq; done

countkmer21.sh GENOME10X.fna       ; mv GENOME10X.fna.21     FASTA-10X-single.21
countkmer21.sh NC_017625.fna       ; mv NC_017625.fna.21     FASTA-single.21   
countkmer21.sh SRR033548.fastq     ; mv SRR033548.fastq.21   FASTQ-single.21
countkmer21.sh SRR033548.fastq.gz  ; mv SRR033548.fastq.21   FASTQGZ-single.21

#cat NC_017625.fna                | countkmer21.sh > SINGLE-stdin.21   # FAILS fasta not supported
#cat NC_017625.fna NC_017625.fna  | countkmer21.sh > DOUBLE-stdin.21   # FAILS fasta not supported
cat SRR033548.fastq | countkmer21.sh > FASTQ-stdin.21 

countkmer21.sh NC_017625.fna GENOME10X.fna SRR033548.fastq SRR033548A.fastq.gz 

mv NC_017625.fna.21      FASTA-multiple.21
mv GENOME10X.fna.21      FASTA-10X-multiple.21
mv SRR033548.fastq.21    FASTQ-multiple.21 
mv SRR033548A.fastq.21   FASTQGZ-multiple.21 

# TEST fastq files big enough to be split:
countkmer21.sh BIG.fastq ; mv BIG.fastq.21 BIG-single.21
cat BIG.fastq | countkmer21.sh > BIG-stdin.21 
countkmer21.sh GENOME10X.fna BIG.fastq 
mv BIG.fastq.21 BIG-multiple.21

# Test empty input
echo -n | countkmer21.sh > empty1.21
echo | countkmer21.sh > empty2.21 


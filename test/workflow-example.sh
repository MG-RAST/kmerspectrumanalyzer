#!/bin/bash
# Download compressed FASTQ from ERA 
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001665/SRR001665_1.fastq.gz > SRR001665_1.fastq.gz 
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001665/SRR001665_2.fastq.gz > SRR001665_2.fastq.gz 

# Stream both mate files to countkmer21.sh 
zcat SRR001665_1.fastq.gz SRR001665_2.fastq.gz | countkmer21.sh > SRR001665-both.21

# Alternatively, countkmer21.sh can be invoked with filenames:
gunzip SRR001665_?.fastq.gz
# This will create SRR001665_1.fastq.21 and 001665_2.fastq.21
countkmer21.sh SRR001665_?.fastq 

# Make visualizations and populate kmer.log

plotkmerspectrum.py SRR001665-both.21  -w png -g 1
plotkmerspectrum.py SRR001665-both.21  -w png -g 5
plotkmerspectrum.py SRR001665-both.21  -w png -g 6

# Fit for genome size 
kmerspectrumanalyzer.py SRR001665-both.21



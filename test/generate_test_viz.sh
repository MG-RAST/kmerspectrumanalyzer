#!/bin/bash

# retrieve 454 DH1 sequencing run, in three parts
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR016/SRR016860/SRR016860.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR016/SRR016859/SRR016859.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR016/SRR016858/SRR016858.fastq.gz
zcat SRR0168??.fastq.gz  > SRR016860A.fastq

# retrieve reference DH1 genome
wget ftp://ftp-trace.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_DH1_uid161951/NC_017625.fna

# generate kmer histogram for sequencing
countkmer21.sh SRR016860A.fastq  
countkmer15.sh SRR016860A.fastq  

# generate kmer histogram for genome
countkmer21.sh NC_017625.fna 

# generates test visualizations
kmerdriver.sh SRR016860A.fastq.21    

# fit 454 data for genome size
plot-kmer-spectrum.py NC_017625.fna.21  SRR016860A.fastq.21 -o NC_017625


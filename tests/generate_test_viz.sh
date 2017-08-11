#!/bin/bash

# retrieve 454 DH1 sequencing run, in three parts
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR016/SRR016860/SRR016860.fastq.gz > SRR016860.fastq.gz
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR016/SRR016859/SRR016859.fastq.gz > SRR016859.fastq.gz
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR016/SRR016858/SRR016858.fastq.gz > SRR016858.fastq.gz
zcat SRR0168??.fastq.gz  > SRR016860A.fastq

# retrieve reference DH1 genome
curl ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/ASSEMBLY_BACTERIA/Escherichia_coli/GCF_000023365/NC_017625.fna > NC_017625.fna

# generate kmer histogram for sequencing
countkmer21.sh SRR016860A.fastq  
countkmer15.sh SRR016860A.fastq  

# generate kmer histogram for genome
countkmer21.sh NC_017625.fna 

# generates test visualizations
kmerdriver.sh SRR016860A.fastq.21    

# fit 454 data for genome size
plotkmerspectrum.py NC_017625.fna.21  SRR016860A.fastq.21 -o NC_017625




# retrieve 454 DH1 sequencing run, in three parts
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR016/SRR016860/SRR016860.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR016/SRR016859/SRR016859.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR016/SRR016858/SRR016858.fastq.gz
gunzip SRR0168??.fastq.gz


wget ftp://ftp-trace.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_DH1_uid161951/NC_017625.fna

# concatenate 
cat SRR016858.fastq SRR016859.fastq SRR016860.fastq > SRR016860A.fastq
countkmer21.sh SRR016860.fastq  # out of memory on micro instance
countkmer15.sh SRR016860.fastq  # out of memory on micro instance

# modified countkmer21.sh  to use smaller initial table allocations -s 40000000
countkmer21.sh SRR016860A.fastq  
countkmer15.sh SRR016860A.fastq  

countkmer21.sh NC_017625.fna 

kmerdriver.sh SRR016860A.fastq.21    # generates test visualizations

plot-kmer-spectrum.py NC_017625.fna  -o NC_017625


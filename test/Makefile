# Makefile to download some data, count some kmers, and produce some visualizations

test: SRR016860A.fastq.21.fit.csv NC_017625.fna.21 SRR016860A.fastq.21
	echo

# retrieve 454 DH1 sequencing run, in three parts
SRR016860.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR016/SRR016860/SRR016860.fastq.gz # 114Mb 
SRR016859.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR016/SRR016859/SRR016859.fastq.gz # 1.6Mb
SRR016858.fastq.gz:
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR016/SRR016858/SRR016858.fastq.gz # 133Mb
SRR016860A.fastq: SRR016860.fastq.gz SRR016858.fastq.gz SRR016859.fastq.gz
	zcat SRR0168??.fastq.gz  > SRR016860A.fastq

# retrieve reference DH1 genome
NC_017625.fna:
	wget ftp://ftp-trace.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_DH1_uid161951/NC_017625.fna

%.fastq.21: %.fastq
	countkmer21.sh $<

%.fna.21: %.fna
	countkmer21.sh $<

%.21.fit.csv: %.21
	-kmerspectrumanalyzer.py $< -n 10 -q

# generates test visualizations
%.21.6.png: %.21
	kmerdriver.sh $^

# fit 454 data for genome size
NC_017625.6.pdf: 
	plot-kmer-spectrum.py NC_017625.fna.21  SRR016860A.fastq.21 -o NC_017625

clean:
	rm *.fit.csv *.fit.detail.csv *.fit.png *.pdf

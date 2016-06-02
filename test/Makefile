# Makefile to download some data, count some kmers, and produce some visualizations

test: SRR016860A.fastq.21.fit.csv NC_017625.fna.21 SRR016860A.fastq.21 NC_017625.fna.6.pdf
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

SRR033548.fastq:
	curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR033/SRR033548/SRR033548.fastq.gz > SRR033548.fastq.gz   # 3Mb tiny fastq-format dataset

SRR001665_1.fastq.gz:
	curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001665/SRR001665_1.fastq.gz > SRR001665_1.fastq.gz   # 320Mb
SRR001665_2.fastq.gz:
	curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001665/SRR001665_2.fastq.gz > SRR001665_2.fastq.gz   # 350Mb

SRR001665_1_X10.fastq: SRR001665_1.fastq
	cat SRR001665_1.fastq SRR001665_1.fastq SRR001665_1.fastq SRR001665_1.fastq SRR001665_1.fastq SRR001665_1.fastq SRR001665_1.fastq SRR001665_1.fastq SRR001665_1.fastq SRR001665_1.fastq > SRR001665_1_X10.fastq

%: %.gz
	cat $^ | gunzip > $@

# retrieve reference DH1 genome
NC_017625.fna:
	wget ftp://ftp.ncbi.nih.gov/genomes/archive/old_refseq/ASSEMBLY_BACTERIA/Escherichia_coli/GCF_000023365/NC_017625.fna  # 4.5Mb

%.fastq.21: %.fastq
	countkmer21.sh $<

%.fna.21: %.fna
	countkmer21.sh $<

%.21.fit.csv: %.21
	-kmerspectrumanalyzer.py $< -n 10 -q

# generates test visualizations
%.21.6.png: %.21
	kmerdriver.sh $^

# display some data
%.6.pdf: %.21
	$(eval STEM=$(subst .21,,$<))
	for i in $$(seq -1 17); do plotkmerspectrum.py $< -o $(STEM) -g $$i; done 

clean:
	rm *.fit.csv *.fit.detail.csv *.fit.png *.pdf


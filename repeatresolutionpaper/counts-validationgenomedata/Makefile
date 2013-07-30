# Makefile for automating the counting and fitting of 21-mers from fastq files from SRA
# The target all counts 21-mers for all the 

FASTQFILES = $(wildcard *.fastq)
HISTFILES = $(subst .fastq,.fastq.21,$(FASTQFILES))
FITFILES  = $(subst .fastq,.fastq.21.fit.csv,$(FASTQFILES))

.PRECIOUS: %.21

all: $(FITFILES)
        echo $(FITFILES)

countall: $(HISTFILES)
        echo $(HISTFILES)

%.fastq.21: %.fastq
#       echo create $@ from $<
        countkmer21.sh $<

%.21.fit.csv : %.21
#       echo create $@ from $<
        -kmerspectrumanalyzer.py $< -n 30 -q

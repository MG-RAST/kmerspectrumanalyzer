#!/bin/bash 
plotkmerspectrum.py ../repeatresolutionpaper/counts-validationgenomedata/SRR000333.fastq.21 -g 25
plotkmerspectrum.py ../repeatresolutionpaper/counts-validationgenomedata/SRR000333.fastq.21 -g 26
cp ../repeatresolutionpaper/counts-validationgenomedata/SRR000333.fastq.21.26.pdf .
cp ../repeatresolutionpaper/counts-validationgenomedata/SRR000333.fastq.21.25.pdf .

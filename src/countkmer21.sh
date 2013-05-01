#!/bin/bash
filename=$1
kmer-tool2  -t fastq -l 21 -f histo -i $1 -o $1.21 

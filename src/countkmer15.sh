#!/bin/bash
filename=$1
k=15
kmer-tool2  -t fastq -l $k -f histo -i $1 -o $1.$k 

#!/bin/bash
filename=$1
k=15
kmer-tool2  -t fastq -l $k -f histo -i $filename -o $filename.${k}

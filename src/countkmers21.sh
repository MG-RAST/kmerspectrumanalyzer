#!/bin/bash
filename=$1
k=21
kmer-tool2  -t fastq -l $k -f dump -i $filename -o $filename.${k}mers

#!/bin/bash
filename=$1
FILETYPE=fastq

if [[ $filename == *".fasta" ]]
then
FILETYPE=fasta
fi
if [[ $filename == *".fna" ]]
then
FILETYPE=fasta
fi
if [[ $filename == *".fa" ]]
then
FILETYPE=fasta
fi

k=21
kmer-tool2  -t $FILETYPE -l $k -f dump -i $filename -o $filename.${k}mers

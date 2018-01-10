#!/bin/bash
filename=$1
k=15
filetype=fastq
USAGE="Usage: countkmer${k}.sh <fastq filename>"
if [ $# -lt 1 ]
    then
    echo "Error: fastq filename is required"
    echo $USAGE 
    exit
elif [ $# -gt 1 ]
    then 
    echo "Error: only one argument is required"
    echo $USAGE 
    exit
elif [[ ! -e $filename ]] 
    then 
    echo "Error: Input filename $filename does not exist."
    echo $USAGE 
    exit
fi 
echo "Counting ${k}mers in $filename, creating $filename.${k} with counts."
kmer-tool2  -t $filetype -l $k -f histo -i $filename -o $filename.${k}

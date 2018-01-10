#!/bin/bash
filename=$1
k=21
# Files are auto-detected, but format is needed in advance
# when streaming data from standard in.  

FILETYPE=fastq   
USAGE1="Usage: countkmer${k}.sh one.fastq
Usage: countkmer${k}.sh one.fasta
       countkmer${k}.sh one.fastq [two.fastq ...]
       countkmer${k}.sh one.fasta [two.fasta ...]
       cat one.fastq two.fastq | countkmer${k}.sh > onetwo.21"


if [ $# -lt 1 ]
    then
    if [ -t 0 ]   # Test STDIN is interactive
        then
        echo "Error: with no arguments, expects FASTQ data from stdin"
        echo "$USAGE1"
        exit
        fi
#        echo "Standard in to default.${k}"
        kmer-tool2  -t $FILETYPE -l $k -f histo -i - -o .default.${k}
        cat .default.${k}
        rm .default.${k}
else  # one or more argument
    for filename in $@ 
        do
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
        if [[ ! -e $filename ]] 
            then 
            echo "Error: Input filename $filename does not exist."
            echo $USAGE 
            exit
#        elif [[ "$filename" =~ ".*.gz" ]]   # doesn't work on some machines?
        elif [[ $( echo $filename | grep .gz$ ) ]] 
                then 
                filestem=${filename/.gz/}
                echo "Counting ${k}mers in $filename, creating $filestem.${k}"
                zcat $filename | kmer-tool2  -t $FILETYPE -l $k -f histo -i - -o $filestem.${k}
        else 
                echo "Counting ${k}mers in $filename, creating $filename.${k}"
                kmer-tool2  -t $FILETYPE -l $k -f histo -i $filename -o $filename.${k}
        fi 
    done
fi

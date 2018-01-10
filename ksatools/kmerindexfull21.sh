#!/bin/bash
stem=$1

if [ ! -e $stem.21.mers ]
then
kmer-tool2  -i $stem -o $stem.21 -l 21 -f dump
cat $stem.21 | sort -n -r -k 2 > $stem.21.mers
else
echo "File $stem.21mers already exists, skipping"
fi

if [ ! -e $stem.21.histhist ] 
then
cut -f 2 $stem.21.mers | sort -n | uniq -c | awk '{print $2 "\t" $1;}' > $stem.21.histhist
else
echo "File $stem.histhist already exists, skipping"
fi 

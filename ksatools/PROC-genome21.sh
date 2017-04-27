#!/bin/bash
stem=$1
out="${stem/.fastq/}"
n=21
mers=$1.$n.mers
scores=$stem.$n.2dscores
twodee=$stem.$n.2d
hi=$stem.hi.fastq

if [ ! -e $mers ]
then
echo running kmerindexfull21.sh   $stem
kmerindexfull21.sh $stem
else
echo "File $mers already exists, skipping"
fi

if [ ! -e $scores ] 
then
echo running fqlookup -1 $stem -i $mers     output $scores
fqlookup.py -1 $stem  -i $mers  > $scores
else
echo File $scores already exists, skipping
fi

if [ ! -e $hi ]
then
echo running fqfilter -1 $stem -i $mers
fqfilter.py -1 $stem  -i $mers   -l -1
else
echo File $scores already exists, skipping
fi

if [ ! -e $twodee ] 
then
echo running join2d.pl
join2d.pl -i 1 -j 2 -s 0 $scores > $twodee
else 
echo File $twodee already exists, skipping
fi


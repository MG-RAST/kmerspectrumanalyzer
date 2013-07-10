#!/bin/bash
target=$1

if [ -e $target ] 
then
plot-kmer-spectrum.py -g 1  $target -w png -o $target
plot-kmer-spectrum.py -g 5  $target -w png -o $target
plot-kmer-spectrum.py -g 6  $target -w png -o $target
plot-kmer-spectrum.py -g 3  $target -w png -o $target
plot-kmer-spectrum.py -g 4  $target -w png -o $target
if [ ! -e $target.10.fit.png ] 
then
kmerspectrumanalyzer.py $target -n 10 -o $target.10
fi
else
echo Can\'t find kmer spectrum $target !
fi
echo "<HTML><HEAD></HEAD><BODY>" > $target.kmers.html
echo "<P><IMG src=$target.1.png>"  >> $target.kmers.html
echo "<P><IMG src=$target.3.png>"  >> $target.kmers.html
echo "<P><IMG src=$target.5.png>"  >> $target.kmers.html
echo "<P><IMG src=$target.6.png>"  >> $target.kmers.html
echo "<P><IMG src=$target.10.fit.png>"  >> $target.kmers.html

echo "<PRE>"   >> $target.kmers.html
cat $target.10.fit.csv >> $target.kmers.html
echo "</PRE>"   >> $target.kmers.html


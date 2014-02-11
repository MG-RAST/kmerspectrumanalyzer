#!/bin/bash
# Run plot-kmer-spectrum on very simple test cases

rm kmers.log

plot-kmer-spectrum.py -g -1  test00.21  # only one count, H = H2 = 1
plot-kmer-spectrum.py -g -1  test01.21
plot-kmer-spectrum.py -g -1  test02.21
plot-kmer-spectrum.py -g -1  test03.21  # uniform distribution  H = H2 = 1000
plot-kmer-spectrum.py -g -1  test04.21
plot-kmer-spectrum.py -g -1  test05.21
plot-kmer-spectrum.py -g -1  test06.21  # all at the same abundance
plot-kmer-spectrum.py -g -1  test07.21  # all at the same abundance 
plot-kmer-spectrum.py -g -1  test08.21  # all at the same abundance
plot-kmer-spectrum.py -g -1  test09.21
plot-kmer-spectrum.py -g -1  test10.21  # rows out of order
plot-kmer-spectrum.py -g -1  test11.21  # zero in the 1st column
plot-kmer-spectrum.py -g -1  test12.21  # zero in the 2nd column
# These should fail to produce any output
plot-kmer-spectrum.py -g -1  testF1.21  # too many zeros
plot-kmer-spectrum.py -g -1  testF2.21  # zero in the 2nd column
plot-kmer-spectrum.py -g -1  testF3.21  # zero in the 1st column
plot-kmer-spectrum.py -g -1  testF4.21  # empty file

if [[ -z $( diff kmers.log kmers-1.log) ]]
then
echo Results match expectation
else
echo kmers.log differs from expectation in kmers-1.log
fi



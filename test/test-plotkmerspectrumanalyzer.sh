#!/bin/bash
rm kmers.log

plot-kmer-spectrum.py -g -1  test0  # only one count, H = H2 = 1
plot-kmer-spectrum.py -g -1  test1
plot-kmer-spectrum.py -g -1  test2
plot-kmer-spectrum.py -g -1  test3  # uniform distribution  H = H2 = 1000
plot-kmer-spectrum.py -g -1  test4
plot-kmer-spectrum.py -g -1  test5
plot-kmer-spectrum.py -g -1  test6  # all at the same abundance
plot-kmer-spectrum.py -g -1  test7  # all at the same abundance 
plot-kmer-spectrum.py -g -1  test8  # all at the same abundance
plot-kmer-spectrum.py -g -1  test9
plot-kmer-spectrum.py -g -1  test10  # rows out of order
plot-kmer-spectrum.py -g -1  testF1  # too many zeros
plot-kmer-spectrum.py -g -1  testF2  # zero in the 2nd column
plot-kmer-spectrum.py -g -1  testF3  # zero in the 1st column
plot-kmer-spectrum.py -g -1  testF4  # empty file
plot-kmer-spectrum.py -g -1  test11  # zero in the 1st column
plot-kmer-spectrum.py -g -1  test12  # zero in the 2nd column

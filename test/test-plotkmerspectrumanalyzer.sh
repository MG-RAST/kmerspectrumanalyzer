#!/bin/bash
# Run plotkmerspectrum on very simple test cases

rm kmers.log

plotkmerspectrum.py -g -1  test00.21  # only one count, H = H2 = 1
plotkmerspectrum.py -g -1  test01.21
plotkmerspectrum.py -g -1  test02.21
plotkmerspectrum.py -g -1  test03.21  # uniform distribution  H = H2 = 1000
plotkmerspectrum.py -g -1  test04.21
plotkmerspectrum.py -g -1  test05.21
plotkmerspectrum.py -g -1  test06.21  # all at the same abundance
plotkmerspectrum.py -g -1  test07.21  # all at the same abundance 
plotkmerspectrum.py -g -1  test08.21  # all at the same abundance
plotkmerspectrum.py -g -1  test09.21
plotkmerspectrum.py -g -1  test10.21  # rows out of order
plotkmerspectrum.py -g -1  test11.21  # zero in the 1st column
plotkmerspectrum.py -g -1  test12.21  # zero in the 2nd column
# These should fail to produce any output
plotkmerspectrum.py -g -1  testF1.21  # too many zeros
plotkmerspectrum.py -g -1  testF2.21  # zero in the 2nd column
plotkmerspectrum.py -g -1  testF3.21  # zero in the 1st column
plotkmerspectrum.py -g -1  testF4.21  # empty file

if [[ -z $( diff kmers.log kmers-1.log) ]]
then
echo Results match expectation
else
echo kmers.log differs from expectation in kmers-1.log
fi



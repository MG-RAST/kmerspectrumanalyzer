#!/bin/bash
# Run plotkmerspectrum on very simple test cases

rm kmers.log

plotkmerspectrum.py test[0-9]?.21  
plotkmerspectrum.py -l testlist1
plotkmerspectrum.py -l testlist2
plotkmerspectrum.py -l testlist3
plotkmerspectrum.py -l testlist4
plotkmerspectrum.py -l testlist5
plotkmerspectrum.py -l testlist6
plotkmerspectrum.py -l testlist7
plotkmerspectrum.py -l testlist8
plotkmerspectrum.py .boguswc*
plotkmerspectrum.py emptyfile 
plotkmerspectrum.py -l mgrlist -t mgm

if [[ -z $( diff kmers.log kmers-2.log) ]]
then
echo Results match expectation
else
echo kmers.log differs from expectation in kmers-2.log
fi 

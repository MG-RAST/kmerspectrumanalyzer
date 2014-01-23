#!/bin/bash
# Run plot-kmer-spectrum on very simple test cases

rm kmers.log

plot-kmer-spectrum.py test[0-9]?.21  
plot-kmer-spectrum.py -l testlist1
plot-kmer-spectrum.py -l testlist2
plot-kmer-spectrum.py -l testlist3
plot-kmer-spectrum.py -l testlist4
plot-kmer-spectrum.py -l testlist5
plot-kmer-spectrum.py -l testlist6
plot-kmer-spectrum.py -l testlist7
plot-kmer-spectrum.py .boguswc*
plot-kmer-spectrum.py emptyfile 
plot-kmer-spectrum.py -l mgrlist -t mgm

#kmerspectrumanalyzer
***

## Purpose
This package contains scripts that summarize, visualize, and 
interpret the kmer spectrum (the histogram of abundances of 
oligonucleotide patterns of fixed length) of short-read 
sequence datasets.  

## Prerequisites
This package depends on numpy, scipy, matplotlib, and 
the University of Maryland's Jellyfish kmer counting library.

*   numpy http://www.numpy.org/
*   scipy http://www.scipy.org/
*   matplotlib http://www.matplotlib.org/
*   jellyfish http://www.cbcb.umd.edu/software/jellyfish/

## License
kmerspectrumanalyzer is under the BSD license; see LICENSE.
Distribution, modification and redistribution, incorporation
into other software, and pretty much everything else is allowed.

## Organization
*   src    -- contains scripts
*   pfge_analysis  -- PFGE gel images (and analysis once generated)
*   repeatresolutionpaper  -- contains data supporting the paper

## Paper
A paper describing interpreting the results of kmercounting has
been published August 2013 in BMC Genomics. 2013 14(1):537
"Rapid quantification of sequence repeats to resolve the size, structure and contents of bacterial genomes."
Williams D, Trimble WL, Shilts M, Meyer F, and Ochman H. PMID: 23924250
The manuscript can be found in 
repeatresolutionpaper/manuscript.

## Authors
*   Will Trimble (Argonne National Laboratory)
*   Travis Harrison (Argonne National Laboratory)
*   David Williams (Yale University)

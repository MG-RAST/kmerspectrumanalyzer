#!/bin/bash

# generates three comparative visualizations from data 
# elsewhere in the repository from lists of files and 
# labels in filelistcb and filelistsz
plotkmerspectrum.py -l filelistcv -g 1 -w png
plotkmerspectrum.py -l filelistsz -g 5 -w png
plotkmerspectrum.py -l filelistsz -g 6 -w png


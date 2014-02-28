#!/bin/bash

# generates all the visualizations

for i in $( seq -1 17 )
do 
plotkmerspectrum.py -l filelistcv -g $i -w png  
done
plotkmerspectrum.py -l filelistcv -g 25 -w png  
plotkmerspectrum.py -l filelistcv -g 26 -w png  


#!/bin/bash

python -c 'import numpy' 
ERR=$?
if [ $ERR -ne 0 ]
then  
echo ERROR: numpy missing
echo See if you can get numpy from http://www.numpy.org/
fi

python -c 'import scipy' 
ERR=$?
if [ $ERR -ne 0 ]
then
echo ERROR: scipy missing
echo See if you can get scipy from http://www.scipy.org/
fi

python -c 'import matplotlib' 
ERR=$?
if [ $ERR -ne 0 ]
then
echo ERROR: matplotlib missing
echo See if you can get matplotlib from http://www.matplotlib.org/
fi

jellyfish -h > /dev/null
ERR=$?
if [ $ERR -ne 0 ]
then
echo ERROR: jellyfish missing
echo See if you can get jellyfish from http://www.cbcb.umd.edu/software/jellyfish/
fi

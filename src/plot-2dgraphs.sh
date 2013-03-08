#!/bin/bash
# script to generate several plots of 2d-histogram data

# Usage:  plot2.sh <matrix.2d>  <xaxislabel> <yaxislabel>
in=$1
in=${in/.2d/}

xaxislabel=$2
yaxislabel=$3

showall.py $in.2d    -o $in.all -x $xaxislabel -y $yaxislabel
showall.py $in.2d    -o $in.100 -x $xaxislabel -y $yaxislabel -m 0 -n 0
showall.py $in.2d    -o $in.010 -x $xaxislabel -y $yaxislabel -m 10 -n 10

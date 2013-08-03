#!/bin/bash

sudo apt-get update
sudo apt-get install -y git 
git clone http://github.com/wltrimbl/kmerspectrumanalyzer.git   # The development versison
# git clone http://github.com/MG-RAST/kmerspectrumanalyzer.git 

echo "export PATH=$PATH:~/kmerspectrumanalyzer/src" >> ~/.bash_profile 
source ~/.bash_profile 

# The ubuntu jellyfish package is jellyfish 1.1.6 , probably only changes are bugfixes.  
# matplotlib needed for visualizations
# scipy needed for fitting
sudo apt-get install -y jellyfish python-matplotlib python-scipy



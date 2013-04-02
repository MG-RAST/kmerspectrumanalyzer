#!/bin/bash

sudo apt-get install git -f
git clone http://github.com/MG-RAST/kmerspectrumanalyzer.git 
echo "export PATH=$PATH:~/kmerspectrumanalyzer/src" >> ~/.bash_profile 
export PATH=$PATH:~/kmerspectrumanalyzer/src

# This fetches jellyfish 1.1.6 , probably only changes are bugfixes.  
sudo apt-get install -f jellyfish  


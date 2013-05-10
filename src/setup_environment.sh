#!/bin/bash
sudo apt-get update
sudo apt-get install -y git 
git clone http://github.com/wltrimbl/kmerspectrumanalyzer.git   # The development versison
echo "export PATH=$PATH:~/kmerspectrumanalyzer/src" >> ~/.bash_profile 
source ./.bash_profile 

# This fetches jellyfish 1.1.6 , probably only changes are bugfixes.  
sudo apt-get install -y jellyfish  
# git clone http://github.com/MG-RAST/kmerspectrumanalyzer.git 

# Needed for visualizations
sudo apt-get install -y python-matplotlib 

# Needed for fitting
sudo apt-get install -y python-scipy


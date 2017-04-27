#!/bin/bash
# This script should be sufficient to install the dependencies on a blank EC2 node

sudo apt-get update && sudo apt-get install -y git make curl
git clone http://github.com/wltrimbl/kmerspectrumanalyzer.git   # The development versison
# git clone http://github.com/MG-RAST/kmerspectrumanalyzer.git

echo export PATH=\$PATH:~/kmerspectrumanalyzer/src >> ~/.bash_profile
source ~/.bash_profile

# The ubuntu jellyfish package is jellyfish 1.1.10
sudo apt-get install -y jellyfish
# matplotlib needed for visualizations
# scipy needed for fitting
sudo pat-get install -y python-matplotlib python-scipy



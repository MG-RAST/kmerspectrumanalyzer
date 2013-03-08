#!/bin/bash

sudo apt-get install git -f
git clone https://github.com/wltrimbl/counting-ecoli.git
echo "export PATH=$PATH:~/counting-ecoli/src" >> ~/.bash_profile 
export PATH=$PATH:~/counting-ecoli/src

#  sudo apt-get install bioperl -f    # NOT NEEDED 

# This fetches jellyfish 1.1.6 , probably only changes are bugfixes.  
sudo apt-get install -f jellyfish  


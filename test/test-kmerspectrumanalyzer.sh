#!/bin/bash

# This runs fits on unrealisitc test data

rm fak-0.fit.*

kmerspectrumanalyzer.py fak-0 -n 1
kmerspectrumanalyzer.py fak-1 -n 1
kmerspectrumanalyzer.py fak-2 -n 1
kmerspectrumanalyzer.py fak-3 -n 1
kmerspectrumanalyzer.py ../repeatresolutionpaper/counts-validationgenomedata/SRR006330.fastq.21

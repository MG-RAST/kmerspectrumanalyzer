#!/usr/bin/env python

import sys

from setuptools import setup

setup(name='ksatools',
      version='1.0',
      description='Kmer spectrum analyzer',
      author='W Trimble',
      author_email='trimble@anl.gov',
      url='https://github.com/wltrimbl/kmerspectrumanalyzer',
      packages=['ksatools'],
      scripts=['scripts//kmerspectrumanalyzer.py', 'scripts/plotkmerspectrum.py', 
               'scripts/countkmer15.sh', 'scripts/countkmer21.sh', 'scripts/countkmers15.sh', 
               'scripts/countkmers21.sh', 'scripts/kmer-tool2', 'scripts/rarify.py' ],
      install_requires=  ['numpy >= 1.6', 'matplotlib >= 1.3', 'scipy >= 0.13']
     )


if sys.version_info < (2, 7):
      sys.stderr.write('ERROR: ksatools requires Python 2.7 or later')


#!/usr/bin/env python

from distutils.core import setup

setup(name='ksatools',
      version='1.0',
      description='Kmer spectrum analyzer',
      author='W Trimble',
      author_email='trimble@anl.gov',
      url='https://github.com/wltrimbl/kmerspectrumanalyzer',
      packages=['ksatools'],
     )

if sys.version_info < (2, 7):
      sys.stderr.write("ERROR: ksatools requires Python 2.7 " +


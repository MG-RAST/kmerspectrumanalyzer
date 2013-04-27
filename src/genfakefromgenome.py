#!/usr/bin/env python

import sys, os
from scipy import stats
import numpy as np
from optparse import OptionParser

def count(x):
   '''Adds list of elements to a dictionary counting the occurrences of each'''
   hash = {}
   for element in x:
      try: 
          hash[element] += 1
      except KeyError:
          hash[element] = 1
   return hash

def nbinompdf(xvalues, poissonlambda, alpha):
    '''Negative binomial with lambda, alpha parameterization.'''
    return np.exp(stats.nbinom.logpmf( xvalues, 1/alpha, 1 / ( 1 + poissonlambda * alpha )  ) )

def nbinomrvs(poissonlambda, alpha, size=1):
    '''Negative binomial with lambda, alpha parameterization.'''
    return       stats.nbinom.rvs( 1/alpha, 1 / ( 1 + poissonlambda * alpha )  , size=size) 

if __name__ == '__main__':
    usage  = "usage: %prog -i <input sequence file> <coverage>"
    parser = OptionParser(usage)
    parser.add_option("-i", "--input",  dest="infile", default=None, help="Input genome spectrum.")
    parser.add_option("-s", "--shape",  dest="shape", default=.04, help="Shape parameter.")
    
    (opts, args) = parser.parse_args()
    coverage = float(args[0])
    filename = opts.infile
    if not (filename and os.path.isfile(filename) ):
        parser.error("Missing input file" )
    sys.stderr.write("Generating fake data from %s\n"%filename)
    shap = float(opts.shape)
    m = []
    for l in open(opts.infile):  
        l = l.rstrip()
        h = l.split("\t")
        try :
            n   = int(h[0]) 
        except: 
            h = l.split(" ")
            n   = int(h[0]) 
        a_n = int(h[1]) 
        m.extend(nbinomrvs( coverage * n   , shap/n , size=a_n ) )
    h = count(m)
    
    for i in sorted(h.keys()) :
        print i, "\t", h[i]
  
  

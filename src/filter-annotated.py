#!/usr/bin/env python

import sys, os, random, re
from optparse import OptionParser
from Bio import SeqIO

def gccontent(sq):
  l = len(sq)
  gc_count=0
  for ch in sq:
     if ch in ['G', 'C', 'g', 'c']:
       gc_count += 1
  try:
    r = float(float(gc_count) / l)
  except:
    r=0
  return r

if __name__ == '__main__':
  usage  = "usage: %prog -i <file1> -o <rejectfile> \n filters interleaved fastq on regex in header, outputs to std out"
  parser = OptionParser(usage)
  parser.add_option("-i", "--in",  dest="inp", default=None, help="input file")
  parser.add_option("-o", "--out",  dest="outp", default=None, help="reject fastq file")
  parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=True, help="Verbose [default off]")
  
  (opts, args) = parser.parse_args()
  if not (opts.inp and os.path.isfile(opts.inp) ):
    parser.error("Missing input file" )

  in_fh  = open(opts.inp)
  if(opts.outp):
    out_fh= open(opts.outp,"w")

  sys.stderr.write("Processing %s ... "%(opts.inp, ))
  pattern=re.compile("med21mer=(\d*)") 
  records=SeqIO.parse(in_fh, "fastq")
  for seq_record1 in records:
    seq_record2=records.next()
    m1=int(pattern.search(seq_record1.description).group(1))
    m2=int(pattern.search(seq_record2.description).group(1))
    g1=gccontent(str(seq_record1.seq))
    g2=gccontent(str(seq_record2.seq))
#    print "%d\t%d\t%.2f\t%.2f"%(m1, m2,g1,g2)
  if 1:
    if m1 > 75 and m2 > 75 :
          SeqIO.write([seq_record1, seq_record2], sys.stdout , "fastq")
    else:
       if opts.outp:
          SeqIO.write([seq_record1, seq_record2], out_fh , "fastq")
  if opts.verbose: sys.stderr.write("Done. \n")

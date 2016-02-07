#!/usr/bin/env python

import sys, os
from optparse import OptionParser
from Bio import SeqIO
from string import maketrans

CHOMPSIZE = 100

def lesserkmer(s):
    '''returns the lesser of a kmer and its reverse complement'''
    t = revc(s)
    if t < s:
        return t
    else:
        return s

def revc(s):
    '''returns reverse complement of a sequence'''
    intab = "AaCcGgTt"
    outtab = "TtGgCcAa"
    trantab = maketrans(intab, outtab)
    t = s.translate(trantab)[::-1]
    return t

def gccontent(sq):
    '''returns float gc content of sequence'''
    length = len(sq)
    gc_count = 0
    for ch in sq:
        if ch in ['G', 'C', 'g', 'c']:
            gc_count += 1
    try:
        r = float(gc_count) / length
    except:
        r = 0
    return r

def read_index(filename):
    gian = {}
    sys.stderr.write("Processing table %s  ...\n"%(filename,))
    in_idx  = open(filename)
    for l in in_idx:
        if l[0] != "#":
            s = l.rstrip().split()
            gian[s[0]] = int(s[1])
    return gian

def kmerabundance(seq, index):
    '''looks up kmer abundance of each kmer in sequence, returns summary statistics; index is the kmer hash'''
    a = []
    for i in range(0, len(seq) - k):
        word = seq[i:i+k]
        w = lesserkmer(str(word).upper())
        if w.find("N") == -1:
            try:
                a.append(index[w])
            except KeyError:
                a.append(0)
    a.sort()
    try:
        median = a[len(a) / 2]
    except IndexError:
        median = 0
    try:
        minimum = a[0]
    except IndexError:
        minimum = 0
    try:
        maximum = a[-1]
    except IndexError:
        maximum = 0
    try:
        average = float(sum(a)) / len(a)
    except IndexError:
        average = 0
    return  (minimum, median, maximum, average)

if __name__ == '__main__':
    usage = "usage: \n"
    parser = OptionParser(usage)
    parser.add_option("-1", "--one", dest="one", default=None, help="Input file 1")
    parser.add_option("-2", "--two", dest="two", default=None, help="Input file 2 (interleaved if absent)")
    parser.add_option("-i", "--index", dest="index", default=None, help="input index(es), comma delimited ")
    parser.add_option("-t", "--type", dest="typ", default="fastq", help="input datatype (fastq, fasta)")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=True, help="Verbose [default off]")
    parser.add_option("-c", "--contig", dest="contig", action="store_true", default=False, help="Process contig")
    (opts, args) = parser.parse_args()
    typ = opts.typ
    if not opts.typ and opts.one[-1] == "a" or opts.one[-1] == "A" or opts.contig:
        typ = "fasta"
    if not (opts.one):
        parser.error("Missing input filename")
    if not os.path.isfile(opts.one):
        parser.error("Missing input file %s"% opts.one)
    if (opts.two and os.path.isfile(opts.two)):
        in_two = open(opts.two)
    in_one = open(opts.one)
    sys.stderr.write("Reading index...\n")
    indexlist = opts.index.split(",")
    indexes = []
    for i in range(len(indexlist)):
        sys.stderr.write("adding index %d" % i)
        giant = read_index(indexlist[i])    # this is the kmer hash
        indexes.append(giant)
    k = len(next(indexes[0].iterkeys()))
    sys.stderr.write("Done slurping... set k = %d\n" % (k,))
#   Setup paired-read input
    sys.stderr.write("Looping data: \n")
    records1 = SeqIO.parse(in_one, typ)
    if(opts.two):
        records2 = SeqIO.parse(in_two, typ)
    else:
        records2 = records1
    n = 0
    if not opts.contig: 
        for seq_record1 in records1:
            n += 1
            seq_record2 = records2.next()
            seq1 = str(seq_record1.seq)
            seq2 = str(seq_record2.seq)
            decoration = ""
            if (seq1.find("N") == -1 and seq2.find("N") == -1):
                for indexno in range(len(indexlist)):
                    (min1, med1, max1, avg1) = kmerabundance(str(seq1) + "N" + str(seq2), indexes[indexno])
                    gc = gccontent(str(seq1)+"N"+str(seq2))
    #                sys.stdout.write("%d\t%d\n%d\t%d\n" % (min1, med1, max1, avg1))
                    decoration = decoration + " %d %d %d %.1f %.3f" % (min1, med1, max1, avg1, gc)
                seq_record1.description = "%s %s" % (seq_record1.description, decoration)
                seq_record2.description = "%s" % (seq_record2.description,) # decoration)
    
                SeqIO.write([seq_record1, seq_record2], sys.stdout, typ)
        if opts.verbose: sys.stderr.write("Done. \n")
    else:   # decorate contigs 
        for seq_record in records1:
            header = seq_record.description
            for i in range(len(seq_record.seq) - CHOMPSIZE):
                for indexno in range(len(indexlist)):
                    if indexno == 0:
                        print "Scores.%s" % header,
                    seq = seq_record.seq[i: i+CHOMPSIZE] 
                    (min1, med1, max1, avg1) = kmerabundance(seq, indexes[indexno]) 
                    print "%d %d %d %.1f" % (min1, med1, max1, avg1),
                print




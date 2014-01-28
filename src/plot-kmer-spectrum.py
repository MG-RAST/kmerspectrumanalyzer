#!/usr/bin/env python
'''Tool to generate graphs of kmer spectra'''

import sys, os
import numpy as np
import matplotlib as mpl
from optparse import OptionParser

from ksatools import getcolor, cleanlabel, getmgrkmerspectrum, calccumsum, printstats, loadfile, makegraphs

def main(filename, opt=6, label=None, n=0):
    '''loads file and invokes makegraphs and printstats.
    Appends graphics from each file onto the figure.
    opt is a symbol for the graph type;
    n is the serial number of successful traces.'''
    logfh = open(opts.logfile, "a")
    if opts.filetype.upper() == "MGM":
        spectrum = getmgrkmerspectrum(filename, mgrkey=MGRKEY)
    elif opts.filetype == "file":
        spectrum = loadfile(filename)
    else:
        raise ValueError("%s is invalid type (valid types are mgm and file)"
            % opts.filetype)
    if spectrum == None:   # Abort this trace--but try to graph the others
        return n
    if label == None:
        label = filename
    if spectrum.shape[1] > 0:
        spectrum = spectrum[np.lexsort((spectrum[:, 1], spectrum[:, 0]))]
        sys.stderr.write("Making graphs for %s\n" % filename)
        try:
            makegraphs(spectrum, filename, opt, label, n=n, dump=opts.dump)
            sys.stderr.write("Printing stats in logfile %s %d\n" %
                (opts.logfile, n))
            printstats(spectrum, filename, filehandle=logfh, n=n)
            printstats(spectrum, filename, filehandle=sys.stdout, n=n)
            n += 1
        except ValueError:   # This catches no data or defective data
            sys.stderr.write("Error printing stats for %s\n" % filename)
            print "Unexpected error:", sys.exc_info()[0]
    else:
        sys.stderr.write("Error with dataset %s\n" % filename)
    return n

if __name__ == '__main__':
    usage = '''usage: plot-kmer-spectrum.py [options] <datafile> [<datafile2> <datafile3>...]
       plot-kmer-spectrum.py [options] -l <list of targets, labels> '''
    parser = OptionParser(usage)
    parser.add_option("-d", "--dump", dest="dump", action="store_true",
         default=False, help="dump table with outputs ")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
         default=False, help="verbose")
    parser.add_option("-o", "--outfile", dest="outfile", action="store",
         default=None, help="dump table with outputs ")
    parser.add_option("-g", "--graph", dest="option", action="store", type="int",
         default="6", help="Graph number ")
    parser.add_option("-i", "--interactive", dest="interactive", action="store_true",
         default=False, help="interactive mode--draw window")
    parser.add_option("-l", "--list", dest="filelist",
         default=None, help="file containing list of targets and labels")
    parser.add_option("-t", "--type", dest="filetype",
         default="file", help="type for file list (file,mgm)")
    parser.add_option("-w", "--writetype", dest="writetype",
         default="pdf", help="file type for output (pdf,png)")
    parser.add_option("-a", "--appendlogfile", dest="logfile",
         default="kmers.log", help="logfile for summary statistics")

    (opts, args) = parser.parse_args()
    graphtype = opts.option
    writetype = opts.writetype
    if len(args) == 0 and not opts.filelist:
        print "Missing input file argument!"
        sys.exit(usage)
    assert writetype == "png" or writetype == "pdf"

    if opts.outfile:
        imagefilename = "%s.%d.%s" % (opts.outfile, graphtype, writetype)
    elif opts.filelist:
        imagefilename = "%s.%d.%s" % (opts.filelist, graphtype, writetype)
    else:
        imagefilename = "%s.%d.%s" % (args[0], graphtype, writetype)
        sys.stderr.write("Warning, using default filename %s\n" % (imagefilename,))
    # only invoke interactive backend if requested with -i
    # this stabilizes behavior on non-interactive terminals
    if not opts.interactive:
        mpl.use("Agg")
    else:
        mpl.use('TkAgg')
    import matplotlib.pyplot as plt
    if opts.filetype == "mgm":
        try:
            MGRKEY = os.environ["MGRKEY"]
        except KeyError:
            MGRKEY = ""

    graphcount = 0
    # Loop over input identifiers, and skip if main()
    # fails to produce some traces
    if opts.filelist:
        assert os.path.isfile(opts.filelist), "File %s does not exist" % opts.filelist
        IN_FILE = open(opts.filelist, "r")
        for line in IN_FILE:
            if line[0] != "#":
                a = line.strip().split("\t")
                if len(a[0]) > 0:
                    if len(a) == 1:
                        a.append(a[0])
                    sys.stderr.write("%s  %s \n" % (a[0], a[1]))
                    graphcount = main(a[0], graphtype, label=a[1], n=graphcount)
    else:
        for f in args:
            filen = f
            graphcount = main(filen, graphtype, n=graphcount)
    # don't continue if all inputs fail
    assert graphcount > 0, "ERROR: unable to find any data to graph!"
    if graphtype != -1:
        sys.stderr.write("Writing graph into file %s\n" % (imagefilename))
        plt.savefig(imagefilename)
    if opts.interactive:
        plt.show()
    else:
        sys.stderr.write("Use -i to open widow with graph\n")
    if graphcount == 0:
        sys.stderr.write("ERROR:  no data found!\n")

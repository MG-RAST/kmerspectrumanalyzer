#!/usr/bin/env python
'''Tool to generate graphs of kmer spectra'''

import sys, os
import numpy as np
import matplotlib as mpl
from optparse import OptionParser

from ksatools import getcolor, cleanlabel, getmgrkmerspectrum, sortbycp, calccumsum, printstats, loadfile

def makegraphs(spectrum, filename, option=6, label=None, n=0):
    '''Draw graphs, one at a time, and add them to the current plot'''
    (cn, c1, yd, yo, zd, zo, y) = calccumsum(spectrum)
    if label == None:
        tracelabel = cleanlabel(filename)
    else:
        tracelabel = cleanlabel(label)
    # sorted by abundance/coverage
    b = np.flipud(np.sort(spectrum.view('float,float'), order=['f0'],
           axis=0).view(np.float))
    # sorted by size (for contigs)
    c = np.flipud(np.sort(spectrum.view('float,float'), order=['f1'],
           axis=0).view(np.float))
    # sorted by abundance-size product (explained)
    d = sortbycp(spectrum)
    (b_cn, b_c1, b_yd, b_yo, b_zd, b_zo, b_y) = calccumsum(b) # abundance
    (c_cn, c_c1, c_yd, c_yo, c_zd, c_zo, c_y) = calccumsum(c) # size
    (d_cn, d_c1, d_yd, d_yo, d_zd, d_zo, d_y) = calccumsum(d) # explained
    No = b_zo.max()
    Nd = b_zd.max()
    x = np.arange(len(b[:, 0]))                  # rank
    color = getcolor(n)
    if option == 0:
        pA = plt.loglog(b_cn, b_c1, "-", color=color, label=tracelabel)
        pA = plt.loglog(b_cn, b_c1, ".", color=color)
        plt.xlabel("kmer abundance")
        plt.ylabel("number of kmers")
        plt.legend(loc="upper right")
    if option == 0 or option == -1:
        if opts.dump:
            c = np.hstack((cn.reshape((len(cn), 1)),
                (c1.reshape((len(cn), 1)))))
            sys.stderr.write("saving output table in %s.0.plot.csv\n" %
                filename)
            np.savetxt("%s.0.plot.csv" % filename, c, fmt=['%d', '%d'],
                delimiter="\t")
    elif option == 1:
        pA = plt.loglog(cn, cn * c1, "-", color=color, label=tracelabel)
        pA = plt.loglog(cn, cn * c1, '.', color=color)
        plt.xlabel("kmer abundance")
        plt.ylabel("kmers observed")
        plt.legend(loc="upper right")
        plt.grid(1)
        if opts.dump:
            c = np.hstack((cn.reshape((len(cn), 1)),
                ((cn * c1).reshape((len(cn), 1)))))
            sys.stderr.write("saving output table in %s.1.plot.csv\n" %
                filename)
            np.savetxt("%s.1.plot.csv" % filename, c, fmt=['%d', '%d'],
                delimiter="\t")
    elif option == 2:
        pA = plt.loglog(b_zo, b_cn, color=color, label=tracelabel)
        plt.xlabel("cumulative kmers observed")
        plt.ylabel("kmer abundance")
        plt.legend(loc="lower left")
        plt.grid(1)
    elif option == 3:
        pA = plt.semilogy(b_zo / No, b_cn, color=color, label=tracelabel)
        pA = plt.semilogy(b_zo / No, b_cn, '.', color=color)
        plt.xlabel("fraction of observed kmers")
        plt.ylabel("kmer abundance ")
        plt.grid(1)
        plt.legend(loc="lower left")
    elif option == 4: # Fraction of distinct kmers vs abundance  NOT RECOMMENDED
        pA = plt.semilogy(b_zd / Nd, b_cn, color=color, label=tracelabel)
        pA = plt.semilogy(b_zd / Nd, b_cn, '.', color=color)
        plt.xlabel("fraction of distinct kmers")
        plt.ylabel("kmer abundance")
        plt.legend(loc="upper right")
        plt.grid(1)
    elif option == 5:
        pA = plt.semilogx(yd, zo / No, '-', color=color)
        pA = plt.semilogx(yd, zo / No, '.', color=color, label=tracelabel)
        plt.xlabel("kmer rank")
        plt.ylabel("fraction of observed kmers")
        plt.xlim((1, 10**10))
        plt.ylim(0, 1)
        plt.grid(1)
        plt.legend(loc="lower left")
    elif option == 6:
        pA = plt.loglog(b_zd, b_cn, '-', color=color, label=tracelabel)
        pA = plt.loglog(b_zd, b_cn, '.', color=color)
        plt.xlabel("kmer rank")
        plt.ylabel("kmer abundance")
        plt.xlim((1, 10**10))
        plt.ylim(1, 10**7)
        plt.grid(1)
        plt.legend(loc="lower left")
        if opts.dump:
            c = np.hstack((yd.reshape((len(yd), 1)), cn.reshape((len(cn), 1))))
            sys.stderr.write("saving output table in %s.6.plot.csv\n" % filename)
            np.savetxt("%s.6.plot.csv" % filename, c, fmt=['%d', '%d'], delimiter="\t")
    elif option == 7:
        pA = plt.plot(x, c_zd, '.-', color=color, label=tracelabel)
        plt.xlabel("contig size rank")
        plt.ylabel("cuml contig size")
        plt.grid(1)
        plt.legend(loc="upper right")
    elif option == 8:
        pA = plt.plot(x, c_zo / No, '.-', color=color, label=tracelabel)
        plt.xlabel("contig size rank ")
        plt.ylabel("frac data explained ")
        plt.grid(1)
        plt.legend(loc="upper right")
    elif option == 9:
        pA = plt.plot(x, d_zo, '-', color=color, label=tracelabel)
        plt.xlabel("contig explain rank ")
        plt.ylabel("data explained ")
        plt.grid(1)
        plt.legend(loc="upper right")
    elif option == 10:
        pA = plt.plot(x, b_yo / No, '.-', color=color, label=tracelabel)
        plt.xlabel("contig cov rank ")
        plt.ylabel("frac data explained ")
        plt.grid(1)
        plt.legend(loc="upper right")
    elif option == 11:
        pA = plt.plot(x, b_yo, '.-', color=color, label=tracelabel)
        plt.xlabel("contig cov rank ")
        plt.ylabel("data explained (bogo bp) ")
        plt.grid(1)
        plt.legend(loc="upper right")
    elif option == 12:
        pA = plt.plot(c_zd, c_zo, '.-', color=color, label=tracelabel)
        plt.xlabel("cumulative contig size")
        plt.ylabel("data explained (bogo bp) ")
        plt.grid(1)
        plt.legend(loc="upper right")
    elif option == 13:
        pA = plt.plot(x, b_cn, '.-', color=color, label=tracelabel)
        plt.xlabel("contig cov rank ")
        plt.ylabel("kmer abundance")
        plt.grid(1)
        plt.legend(loc="upper right")
    elif option == 14:
        pA = plt.plot(x, d_cn * d_c1, '.-', color=color, label=tracelabel)
        plt.xlabel("contig expl rank ")
        plt.ylabel("data explained (bogo bp) ")
        plt.grid(1)
        plt.legend(loc="upper right")
    elif option == 15:
        pA = plt.plot(x, c_c1, '.-', color=color, label=tracelabel)
        plt.xlabel("contig size rank")
        plt.ylabel("contig size")
        plt.grid(1)
        plt.legend(loc="upper right")
    elif option == 16:
        pA = plt.plot(x, c_yd, '.-', color=color, label=tracelabel)
        plt.xlabel("contig expl rank ")
        plt.ylabel("data explained (bogo bp) ")
        plt.grid(1)
        plt.legend(loc="upper right")
    elif option == 17:
        pA = plt.plot(x, d_cn * d_c1, '.-', color=color, label=tracelabel)
        plt.xlabel("contig expl rank ")
        plt.ylabel("data explained (bogo bp) ")
        plt.grid(1)
        plt.legend(loc="upper right")

def main(filename, opt=6, label=None, n=0):
    '''loads file and invokes makegraphs and printstats.
    Appends graphics from each file onto the figure.
    opt is a symbol for the graph type;
    n is the serial number of successful traces.'''
    logfh = open(opts.logfile, "a")
    if opts.filetype.upper() == "MGM":
        spectrum = getmgrkmerspectrum(filename, MGRKEY=MGRKEY)
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
            makegraphs(spectrum, filename, opt, label, n=n)
            sys.stderr.write("Printing stats in logfile %s %d\n" %
                (opts.logfile, n))
            printstats(spectrum, filename, filehandle=logfh, n=n)
            printstats(spectrum, filename, filehandle=sys.stdout, n=n)
            n += 1
        except:
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

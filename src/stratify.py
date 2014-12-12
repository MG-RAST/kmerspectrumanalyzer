#!/usr/bin/env python
'''Tool to generate bar graphs of data-depth fractions and
basepair-by-depth distributions from kmer spectra'''

import sys, os
import numpy as np
import matplotlib as mpl
from optparse import OptionParser

import ksatools

def plotstratifiedfracs(labels, spectra):
    '''Given list of labels and list of spectra, splits spectra
    up and produces bar plot of the stratified spectra's usage 
    fractions'''
    colors = ["#CCFFCC", "#99FF99", "#55FF55", "#33CC33", "#009933", "#003300"]
    bands = [] ; fracs = []; sizes = [] ; fracc = []
    plt.grid(linestyle="-", zorder=-10)
    for i in range(len(labels)):
        label = labels[i]
        spectrum = spectra[i]
        sys.stderr.write("Stratifying "+label+"...\n")
        BANDS = [1, 3, 30, 300, 3000, 30000, 30000000]
        band, frac, size = ksatools.stratify(spectrum, bands=BANDS)
        bands.append(band)
        fracs.append(frac)
        sizes.append(size)
    for l in range(len(labels)):
        for i in range(len(bands[0])-1):
            if l == 0:
                plt.barh(l, (fracs[l][i]-fracs[l][i+1]), left=(fracs[l][i+1]), color=colors[i], label=str(bands[0][i])+"-"+str(bands[0][i+1]), alpha=1.0, zorder=0)
            else:
                plt.barh(l, (fracs[l][i]-fracs[l][i+1]), left=(fracs[l][i+1]), color=colors[i], alpha=1.0, zorder=0)
    pos = np.arange(len(labels)) + 0.5
    plt.yticks(pos, labels)
    plt.xlabel("Cumulative data fraction")
    plt.tight_layout()
    if not opts.suppresslegend:
        plt.legend(loc="upper left")
    plt.show()

def plotstratifiedsizes(labels, spectra):
    '''Given list of labels and spectra, produces stacked bar graphs
    on a log scale of the cumulative number of basepairs above or
    equal to each depth boundary.'''
    colors = ["#CCFFCC", "#99FF99", "#55FF55", "#33CC33", "#009933", "#003300"]
    plt.grid(linestyle="-", zorder=-10)
    bands = [] ; fracs = []; sizes = [] ; fracc = []
    BANDS = [1, 3, 30, 300, 3000, 30000, 30000000]
    for i in range(len(labels)):
        label = labels[i]
        spectrum = spectra[i]
        sys.stderr.write("Stratifying "+label+"...\n")
        band, frac, size = ksatools.stratify(spectrum, bands=BANDS)
        bands.append(band)
        fracs.append(frac)
        sizes.append(size)
    for l in range(len(labels)):
        for i in range(len(bands[0])-1):
            sizec = np.array(sizes[l])
            if l == 1:
                plt.barh(l, (sizec[i+1]-sizec[i]), left=(sizec[i]), color=colors[i], label=str(bands[0][i])+"-"+str(bands[0][i+1]), log=True)
            else:
                plt.barh(l, (sizec[i+1]-sizec[i]), left=(sizec[i]), color=colors[i], log=True)
    pos = np.arange(len(labels)) + 0.5
    plt.xlim((1, 1E9))
    plt.yticks(pos, labels)
    plt.xlabel("Distinct kmers (basepairs)")
    plt.tight_layout()
    if not opts.suppresslegend:
        plt.legend()
    plt.show()

def summarizestrata(labels, spectra):
    '''Prints one-line table-style summary of cumulative fractions and
    sizes '''
    BANDS = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000]
    for spectrum, label in zip(spectra, labels):
        band, frac, size = ksatools.stratify(spectrum, bands=BANDS)
        bandsbanner = map(str, BANDS)
        print "#name\t" + "\t".join(bandsbanner) + "\t" + "\t".join(bandsbanner)
        print label+"\t" + "\t".join(map(str, frac)) + "\t",
        print "" + "\t".join(map(str, size))
    return

if __name__ == '__main__':
    usage = '''spectra.py [options] -l <list of target files, labels>
'''
    parser = OptionParser(usage)
    parser.add_option(
        "-v", "--verbose", dest="verbose", action="store_true",
        default=False, help="verbose")
    parser.add_option(
        "-o", "--outfile", dest="outfile", action="store",
        default=None, help="output file name")
    parser.add_option(
        "-g", "--graph", dest="option", action="store", type="int",
        default="6", help="graph type 0 = fraction, 1 = basepairs")
    parser.add_option(
        "-i", "--interactive", dest="interactive", action="store_true",
        default=False, help="interactive mode--draw window")
    parser.add_option(
        "-l", "--list", dest="filelist",
        default=None, help="file containing list of targets and labels")
    parser.add_option(
        "-t", "--type", dest="filetype",
        default="file", help="type for file list (file,mgm)")
    parser.add_option(
        "-w", "--writetype", dest="writetype",
        default="pdf", help="file type for output (pdf,png)")
    parser.add_option(
        "-s", "--suppresslegend", dest="suppresslegend", action="store_true",
        default=False, help="supress display of legend")
    parser.add_option(
        "-n", "--name", dest="title",
        default=None, help="Name for graph, graph title")

    (opts, args) = parser.parse_args()
    writetype = opts.writetype
    assert writetype == "png" or writetype == "pdf" or writetype == "eps"
    if opts.outfile:
        imagefilename = opts.outfile
    else:
        imagefilename = opts.filelist + "." + opts.writetype
        sys.stderr.write("Warning, using default filename %s\n" % (imagefilename,))
    # only invoke interactive backend if requested with -i
    # this stabilizes behavior on non-interactive terminals
    if not opts.interactive:
        mpl.use("Agg")
    else:
        mpl.use('TkAgg')
    import matplotlib.pyplot as plt

    graphcount = 0
    # Loop over input identifiers, and skip if main()
    # fails to produce some traces
    spectra = []
    labels = []
    if opts.filelist:
        assert os.path.isfile(opts.filelist), "File %s does not exist" % opts.filelist
        IN_FILE = open(opts.filelist, "r")
        for line in IN_FILE:
            if line[0] != "#":
                a = line.strip().split("\t")
                if len(a[0]) > 0:
                    if len(a) == 1:
                        a.append(a[0])
                    sys.stderr.write("%s\t%s\n" % (a[0], a[1]))
                    spectra.append(ksatools.loadfile(a[0]))
                    labels.append(a[1])
    if opts.option == 0:
        plotstratifiedfracs(labels, spectra)
    if opts.option == 1:
        plotstratifiedsizes(labels, spectra)
    if opts.option == -1:
        summarizestrata(labels, spectra)
        sys.exit()
    sys.stderr.write("Writing graph into file %s\n" % (imagefilename))
    if opts.interactive:
        plt.show()
    else:
        sys.stderr.write("Use -i to open widow with graph\n")
    plt.savefig(imagefilename)

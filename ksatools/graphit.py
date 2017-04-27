#!/usr/bin/env python
'''Tool to generate graphs of kmer spectra'''

COLORLIST = ["b", "g", "r", "c", "y",
             "m", "k", "BlueViolet", "Coral", "Chartreuse",
             "DarkGrey", "DeepPink", "LightPink"]

import sys
import numpy as np
import ksatools
import matplotlib as mpl
from optparse import OptionParser

def getcolor(index, colorlist):
    if colorlist == []:
        colorlist = COLORLIST
    l = index % len(colorlist)
    return colorlist[l]

def plotme(data, graphtype=None, label=None, n=0, opts=None, color=None, style="-", scale=1):
    import matplotlib.pyplot as plt
    # note, calccumsum will raise an exception here if data is invalid
    if color == None:
        color = getcolor(n, colorlist)
    if label == "":
        label = None
    if graphtype == "linear" or graphtype == None:
#        if opts.markers:
#            pA = plt.plot(data[:, 0], data[:, 1], ".", color=color, label=label, linestyle=style)
        pA = plt.plot(s * data[:, 0], data[:, 1], color=color, label=label, linestyle=style)
        legendloc = "upper right"
    if graphtype == "semilogy":
        if opts.dot:
            pA = plt.semilogy(s*data[:, 0], data[:, 1], ".", color=color, label=label, linestyle=style)
        pA = plt.semilogy(s*data[:, 0], data[:, 1], color=color, label=None, linestyle=style, linewidth=opts.thickness)
        legendloc = "upper right"
    if graphtype == "semilogx":
        if opts.dot:
            pA = plt.semilogx(data[:, 0], data[:, 1], ".", color=color, label=label, linestyle=style)
        pA = plt.semilogx(s*data[:, 0], data[:, 1], color=color, label=label, linestyle=style, linewidth=opts.thickness)
        legendloc = "upper right"
    if graphtype == "loglog":
        pA = plt.loglog(s*data[:, 0], data[:, 1], ".", color=color, label=label, linestyle=style)
        pA = plt.loglog(s*data[:, 0], data[:, 1], color=color, label=None, linestyle=style, linewidth=opts.thickness)
        legendloc = "upper right"
    if graphtype == "diff":
        pA = plt.plot(data[1:, 0], np.exp(np.diff(np.log(data[:, 1])))/data[1:, 0], ".", color=color, label=label, linestyle=style)
        pA = plt.plot(data[1:, 0], np.exp(np.diff(np.log(data[:, 1])))/data[1:, 0], color=color, label=Nonte, linestyle=style)
        legendloc = "upper right"
    if not opts.suppress:
        plt.legend()
    if opts.plotlegend is not None:
        plt.gcf().suptitle(opts.plotlegend, fontsize=24, x=0.03)
    if opts.xlim != "":
        x1, x2 = opts.xlim.split(",")
        plt.xlim([float(x1), float(x2)])
    plt.xlabel(opts.xlabel, fontsize=18)
    plt.ylabel(opts.ylabel, fontsize=18)
    plt.grid(1)

if __name__ == '__main__':
    usage = "graphit.py <options> <arguments>"
    parser = OptionParser(usage)
    parser.add_option("-x", "--xlabel", dest="xlabel", action="store",
                      default="x label", help="")
    parser.add_option("-y", "--ylabel", dest="ylabel", action="store",
                      default="y label", help="")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
                      default=False, help="verbose")
    parser.add_option("-o", "--outfile", dest="outfile", action="store",
                      default="test.png", help="dump table with outputs ")
    parser.add_option("-g", "--graphtype", dest="graphtype", action="store",
                      default=None, help="graph type")
    parser.add_option("-i", "--interactive", dest="interactive", action="store_true",
                      default=False, help="interactive mode--draw window")
    parser.add_option("-l", "--list", dest="filelist",
                      default=None, help="file containing list of targets and labels")
    parser.add_option("-t", "--thickness", dest="thickness",
                      default=2, help="line thickness for traces")
    parser.add_option("-w", "--writetype", dest="writetype",
                      default="pdf", help="file type for output (pdf,png)")
    parser.add_option("-p", "--plotlegend", dest="plotlegend",
                      default=None, help="Overall number at top of graph")
    parser.add_option("-s", "--suppresslegend", dest="suppress", action="store_true",
                      default=False, help="supress display of legend")
    parser.add_option("-n", "--name", dest="title",
                      default=None, help="Name for graph, graph title")
    parser.add_option("-c", "--scale", dest="scale",
                      default=False, action="store_true", help="Multiply by col 2")
    parser.add_option("--xlim", dest="xlim",
                      default="", type="str", help="xlimits: comma-separated")
    parser.add_option("--ylim", dest="ylim",
                      default="", type="str", help="ylimits: comma-separated")
    parser.add_option("-d", "--dot", dest="dot",
                      default=False, action="store_true", help="plot dots")

    (OPTS, ARGS) = parser.parse_args()
    SCALE = OPTS.scale
    if not OPTS.interactive:
        mpl.use("Agg")
    else:
        mpl.use('TkAgg')
    import matplotlib.pyplot as plt
    listfile = OPTS.filelist

    IN_FILE = open(listfile, "r").readlines()

    numplots = len(IN_FILE)
    n = 0
    for line in IN_FILE:
        if line[0] != "#":
            a = line.strip().split("\t")
            if len(a[0]) > 0:
                if len(a) == 1:
                    a.append(a[0])
                sys.stderr.write("%s  %s \n" % (a[0], a[1]))
                filename = a[0]
                if len(a) == 3 or len(a) == 4:
                    selectedcolor = a[2]
                else:
                    selectedcolor = getcolor(n, COLORLIST)
                spectrum = ksatools.loadfile(filename)
                if SCALE:
                    s = float(a[1])
                else:
                    s = 1
                if len(a) == 4:
                    selectedcolor = a[2]
                    selectedstyle = a[3]
                    plotme(spectrum, label=a[1], color=selectedcolor, scale=s,
                           style=selectedstyle, graphtype=OPTS.graphtype, opts=OPTS)
                else:
                    plotme(spectrum, label=a[1], color=selectedcolor, scale=s,
                           graphtype=OPTS.graphtype, opts=OPTS)

                    n = n + 1
    if OPTS.suppress == 0:
        plt.legend(loc="upper left")
    else:
        for v in ARGS:
            print(v)
            filename = v
            spectrum = ksatools.loadfile(filename)
            plotme(spectrum, filename, opts=OPTS,
                   color=COLORLIST[n], graphtype=OPTS.graphtype)
            n = n + 1
#        plt.legend(loc="upper left")
    if OPTS.interactive:
        plt.show()
    if OPTS.outfile == "test.png":
        sys.stderr.write("Warning! printing graphs in test.png!\n")
    else:
        sys.stderr.write("Printing graphs in " + OPTS.outfile + "\n")
        plt.savefig(OPTS.outfile)


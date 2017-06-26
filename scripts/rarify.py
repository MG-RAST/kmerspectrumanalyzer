#!/usr/bin/env python
'''Tool to generate computationally-rarefied graphs kmer spectra'''

import numpy as np
import sys, os
import scipy.stats
import ksatools.ksatools 
import matplotlib as mpl
from optparse import OptionParser

from ksatools.rare import fract, rich, calc_resampled_fraction, calc_resampled_richness, plotme 
from ksatools.ksatools import loadfile

if __name__ == "__main__":
    PARSER = OptionParser("rare.py [options] countfilaname [countfilename2...]\n   -- computationally rarefy kmer spectra")
    PARSER.add_option("-i", "--interactive", dest="interactive",
                      action="store_true", default=False,
                      help="interactive mode--draw window")
    PARSER.add_option("-l", "--list", dest="filelist", default=None,
                      help="file containing list of targets and labels")
    PARSER.add_option("-g", "--graphtype", dest="graphtype", default=1,
                      help="graph type 1: shaded 2: non-shaded 3: kmer richness")
    PARSER.add_option("-s", "--suppress", dest="suppresslegend", default=False,
                      action="store_true", help="suppress legend")
    PARSER.add_option("-c", "--colors", dest="colors",
                      help="comma-separated color list")
    PARSER.add_option("-t", "--threshold", dest="threshold", type="int",
                      help="threshold")
    PARSER.add_option("-o", "--output", dest="outfile", default="",
                      help="filename for output")
    PARSER.add_option("-d", "--dump", dest="dump", default=None,
                      action="store_true", help="output table .rare.csv")
    (OPTS, ARGS) = PARSER.parse_args()
    SHADED = int(OPTS.graphtype)
    if len(ARGS) == 0 and not OPTS.filelist:
       sys.stderr.write("Error, requires one or more kmer histogram input filenames.\nrare.py -h lists options\n")
       sys.exit(1)
    n = 0
    if not OPTS.interactive:
        mpl.use("Agg")
    else:
        mpl.use('TkAgg')
    import matplotlib.pyplot as plt
    if OPTS.colors:
        COLORS = OPTS.colors.split(",")
    else:
        COLORS = ["b", "g", "r", "c", "y", "m", "k", "BlueViolet",
                  "Coral", "Chartreuse", "DarkGrey", "DeepPink",
                  "LightPink"]
# construct range of thresholds, calculate threshold fraction curves
# not lightning fast but should be
    listofthresholds = [1, 3.3, 10, 33, 100, 330, 1000, 3300, 10000]
    listofthresholds = 10**np.arange(0, 4.5, 0.5)
    if SHADED == 2 or SHADED == 3:
        listofthresholds = [1]
    else:
        listofthresholds = [1, 3, 10, 30]
    if OPTS.threshold:
        listofthresholds = [OPTS.threshold]
    OUTFILE = OPTS.outfile
    if OUTFILE == "":
        if OPTS.filelist:
            OUTFILE = OPTS.filelist + ".rare." + str(SHADED) + ".png"
        else:
            OUTFILE = "test" + ".rare." + str(SHADED) + ".png"

    if OPTS.filelist:
        listfile = OPTS.filelist
        assert os.path.isfile(listfile), \
            "File {} does not exist".format(listfile)
        IN_FILE = open(listfile, "r").readlines()
        numplots = len(IN_FILE)
        for line in IN_FILE:
            if line[0] != "#":
                a = line.strip().split("\t")
                if len(a[0]) > 0:
                    if len(a) == 1:
                        a.append(a[0])
                    sys.stderr.write("%s  %s \n" % (a[0], a[1]))
                    filename = a[0]
                    if len(a) == 3:
                        selectedcolor = a[2]
                    else:
                        selectedcolor = COLORS[n%len(COLORS)]
                    spectrum = ksatools.loadfile(filename)
                    if spectrum != []:
                        plotme(spectrum, label=a[1], color=selectedcolor,
                               thresholdlist=listofthresholds,
                               numplots=numplots, dump=OPTS.dump, shaded=SHADED)
                        n = n + 1
        if OPTS.suppresslegend == 0:
            plt.legend(loc="upper left")
        plt.savefig(OUTFILE)
    else:
        for v in ARGS:
            print("#", v)
            filename = v
            spectrum = ksatools.loadfile(filename)
            plotme(spectrum, filename, thresholdlist=listofthresholds,
                   color=COLORS[n], dump=OPTS.dump, numplots=len(ARGS), shaded=SHADED)
            n = n + 1
#        plt.legend(loc="upper left")
        sys.stderr.write("Warning! printing graphs in default filename " + OUTFILE + "\n")
        plt.savefig(OUTFILE)

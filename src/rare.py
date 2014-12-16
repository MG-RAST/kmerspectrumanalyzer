#!/usr/bin/env python
'''Tool to generate computationally-rarefied graphs kmer spectra'''

import numpy as np
import matplotlib.pyplot as plt
import sys, os
import scipy.stats
import ksatools
from optparse import OptionParser

def fract(aa, epsilon, threshold):
    '''Evaluates the fraction of theoretically-subsampled spectra
    above a specified threshold.  Dataset abundance is attenuated by
    the factor epsilon.  Returns a float.  aa is a two-column abudnance
    table, epsilon and threshold are floats.'''
    sys.stderr.write("E %f T %f\n" % (epsilon, threshold))
        
    xr = aa[:, 0]
    xn = aa[:, 1]
    NO = np.sum(xn * xr)
    p = 0.0
    smallr = xr * epsilon
    for i in range(len(xr)):
        if smallr[i] > 10 * threshold:
            interim = float(xn[i]*xr[i])
        elif smallr[i] < threshold / 10:
            interim = 0
        else:
            interim = float(xn[i] * xr[i]) * (1 - scipy.stats.binom.cdf(
                threshold + 0.5, xr[i], epsilon)) / (1 - scipy.stats.binom.cdf(
                    0.5, xr[i], epsilon))
        if not np.isnan(interim):
            p += interim
    return p / NO

def calc_resampled_fraction(aa, samplefracs, thresholds):
    '''calculate 2D array of return value of fract by evaluating it
    for each fraction in samplefracs and each threshold in thresholds.
    Returns 2d matrix sith shape = len(samplefracs), len(thresholds)
    aa must be 2d ndarray
    '''
    assert aa.shape[1] == 2
    matrix = np.zeros((len(samplefracs), len(thresholds)))
    for i, frac in enumerate(samplefracs):
        for j, threshold in enumerate(thresholds):
            dummy = fract(aa, frac, threshold)
            matrix[i][j] = dummy
    return matrix

def plotme(b, label, color=None, thresholdlist=None, numplots=4,
     suppress=False, dump=False):
    '''performs calculations and calls graphing routines,
    given spectra'''
# define range of subsamples
    N = np.sum(b[:, 0] * b[:, 1])
    samplefractions = 10**np.arange(2, 11, .5) / N  # CHEAP
    samplefractions = 10**np.arange(2, 11, .1) / N
# Throw away unecessary samples
    samplefractions = np.hstack((samplefractions[samplefractions < 1], 1))
    if thresholdlist == None:
        thresholdlist = [1]

    matrix = calc_resampled_fraction(b, samplefractions, thresholdlist)
    effort = N * samplefractions
    data = np.hstack([np.atleast_2d(effort).T, matrix])
    np.savetxt(sys.stdout, data, fmt="%.3f")
    if dump:
        headertext = "subsetsize\t"+"\t".join(map(str, thresholdlist))
        np.savetxt(label+".rare.csv", data, header=headertext, delimiter="\t")
    pex2 = np.hstack((effort[0], effort, effort[-1]))
    pex = effort
    for i in range(matrix.shape[1]):
        aug2 = np.hstack((0, matrix[:, i], 0))
        aug = matrix[:, i]
#        lab = label + " " + str(thresholdlist[i])
        lab = str(thresholdlist[i]) + "x"
        plt.grid(1)
        if SHADED == 0:
            plt.title(label)
            plt.semilogx(pex, aug, "-o", label=lab)
        elif SHADED == 2:
            lab = label + str(thresholdlist[i]) + "x"
            lab = label
            plt.semilogx(pex, aug, "-", label=lab, color=color)
        elif SHADED == 1:
            plt.subplot(numplots, 1, n + 1)
            plt.semilogx(pex, aug, "-", label=lab, color=color)
            plt.fill(pex2, aug2, "k", alpha=0.2)
            plt.title(label)
        else:
            plt.semilogx(pex, aug, "-", label=lab)
            plt.fill(pex2, aug2, "k", alpha=0.2)
            plt.title(label)
#            label=str(thresholdlist[i]))
#        plt.fill(pex, aug, "k", alpha=0.2)
    plt.ylim((0, 1))
    plt.xlim((1E4, 1E11))
    if SHADED == 0 or n+1 == numplots:
        plt.xlabel("Sequencing effort (bp)")
    else:    # suppress drawing of x-axis labels for all but last plot
        frame1 = plt.gca()
        frame1.axes.get_xaxis().set_ticks([])
    if SHADED == 0 or n == 2 or 1:
        plt.ylabel("Fraction of data")
    plt.tight_layout()
    return()

if __name__ == "__main__":
    PARSER = OptionParser("rare.py -- rarefy kmer spectra")
    PARSER.add_option("-i", "--interactive", dest="interactive",
         action="store_true", default=False,
         help="interactive mode--draw window")
    PARSER.add_option("-l", "--list", dest="filelist", default=None,
         help="file containing list of targets and labels")
    PARSER.add_option("-g", "--graphtype", dest="graphtype", default=1,
         help="graph type 1: shaded 2: non-shaded 3: kmer richness")
    PARSER.add_option("-s", "--suppress", dest="suppresslegend", default=True,
         action="store_true", help="suppress legend")
    PARSER.add_option("-c", "--colors", dest="colors",
         help="comma-separated color list")
    PARSER.add_option("-d", "--dump", dest="dump",
         action="store_true", help="output table .rare.csv")
    (OPTS, ARGS) = PARSER.parse_args()
    SHADED = int(OPTS.graphtype)
    n = 0
    if OPTS.colors:
        COLORS = OPTS.colors.split(",")
    else:
        COLORS = [
                  "b", "g", "r", "c", "y", "m", "k", "BlueViolet",
                  "Coral", "Chartreuse", "DarkGrey", "DeepPink", 
                  "LightPink" ]
# construct range of thresholds, calculate threshold fraction curves
# not lightning fast but should be
    listofthresholds = [1, 3.3, 10, 33, 100, 330, 1000, 3300, 10000]
    listofthresholds = 10**np.arange(0, 4.5, 0.5)
    if SHADED == 2:
        listofthresholds = [1]
    else:
        listofthresholds = [1, 3, 10, 30]
    if OPTS.filelist:
        listfile = OPTS.filelist
        assert os.path.isfile(listfile), "File {} does not exist".format(
             listfile)
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
                    plotme(spectrum, label=a[1], color=selectedcolor,
                        thresholdlist=listofthresholds,
                        numplots=numplots, dump=OPTS.dump)
                    n = n + 1
        if OPTS.suppresslegend == 0:
            plt.legend(loc="upper left")
        plt.savefig(listfile + ".rare."+str(SHADED)+".png")
    else:
        for v in ARGS:
            print v
            filename = v
            spectrum = ksatools.loadfile(filename)
            plotme(spectrum, filename, thresholdlist=listofthresholds,
               color=COLORS[n], dump=OPTS.dump)
            n = n + 1
#        plt.legend(loc="upper left")
        sys.stderr.write("Warning! printing graphs in test."
                         + str(SHADED)+".png!\n")
        plt.savefig("test."+str(SHADED)+".png")

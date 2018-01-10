#!/usr/bin/env python
'''Tool to generate computationally-rarefied graphs kmer spectra'''

import sys
import os
import scipy.stats
from optparse import OptionParser
import numpy as np
import ksatools


def fract(aa, epsilon, threshold):
    '''Evaluates the fraction of theoretically-subsampled spectra
    above a specified threshold.  Dataset abundance is attenuated by
    the factor epsilon.  Returns a float beween 0 and 1.  aa is a
    two-column abudnance table, epsilon and threshold are floats.
    '''
    sys.stderr.write("E %f T %f\n" % (epsilon, threshold))
    xr = aa[:, 0]
    xn = aa[:, 1]
    NO = np.sum(xn * xr)
    p = 0.0
    for i in range(len(xr)):
        # this is the expected number of nonzero categories after hypergeometric sampling
        #        nonzero = (1.-scipy.stats.hypergeom.cdf(0.5, NO, xr[i], epsilon*NO))
        nonzero = (1. - scipy.stats.hypergeom.pmf(0, NO, xr[i], epsilon * NO))
        # For efficiency, don't evaluate if numerator is too small
        # For numerical stability, don't evaluate term if denominator (nonzero) is too small
        # note: second threshold (on nonzero) here creates kinks in the graph, but is important
        if nonzero * xr[i] * xn[i] > 10E-0 and nonzero > 1E-2:
            # and this is the expected number of above-threshold survivors
            gt_thresh = 1. - \
                scipy.stats.hypergeom.cdf(
                    threshold + 0.5, NO, xr[i], epsilon * NO)
            interim = float(xn[i] * xr[i]) * (gt_thresh / nonzero)
            if (not np.isnan(interim)) and (interim > 0):
                p += interim
    return p / NO


def rich(aa, epsilon, threshold):
    sys.stderr.write("richness E %f T %f\n" % (epsilon, threshold))
    xr = aa[:, 0]
    xn = aa[:, 1]
    NO = np.sum(xn * xr)
    interim = 0
    for i in range(len(xr)):
        # this is the expected number of nonzero categories after hypergeometric sampling
        #        nonzero = (1.-scipy.stats.hypergeom.cdf(0.5, NO, xr[i], epsilon*NO))
        nonzero = (1. - scipy.stats.hypergeom.pmf(0, NO, xr[i], epsilon * NO))
        interim += nonzero * xn[i]
    return interim


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


def calc_resampled_richness(aa, samplefracs, thresholds):
    '''calculate 2D array, like calc_resampled_richness, of
    calculated subsampled richness for each fraction in samplefracs
    and each threshold in thresholds.
    Returns 2d matrix sith shape = len(samplefracs), len(thresholds)
    aa must be 2d ndarray
    '''
    assert aa.shape[1] == 2
    matrix = np.zeros((len(samplefracs), len(thresholds)))
    for i, frac in enumerate(samplefracs):
        for j, threshold in enumerate(thresholds):
            dummy = rich(aa, frac, threshold)
            matrix[i][j] = dummy
    return matrix


def plotme(b, label, color=None, thresholdlist=None, numplots=4,
           suppress=False, dump=False, shaded=0, n=1):
    '''Performs calculations and calls graphing routines,
    given spectra
    '''
# define range of subsamples
    import matplotlib.pyplot as plt
    N = np.sum(b[:, 0] * b[:, 1])
    samplefractions = 10**np.arange(2, 11, .5) / N  # CHEAP
    samplefractions = 10**np.arange(2, 11, .1) / N
# Throw away unecessary samples
    samplefractions = np.hstack((samplefractions[samplefractions < 1], 1))
    SHADED = shaded
    if thresholdlist is None:
        thresholdlist = [1]
    if SHADED != 3:
        matrix = calc_resampled_fraction(b, samplefractions, thresholdlist)
    else:
        matrix = calc_resampled_richness(b, samplefractions, thresholdlist)
    effort = N * samplefractions
    data = np.hstack([np.atleast_2d(effort).T, matrix])
#    np.savetxt(sys.stdout, data, fmt="%.3f")   # Numpy can't write to standard out in python3
    headertext = "subsetsize\t" + "\t".join(map(str, thresholdlist))
    with open(label + ".rare.csv", 'wb') as fp:
        np.savetxt(fp, data, header=headertext, delimiter="\t")
    if dump:
        with open(label + ".rare.csv") as f:
            for l in f:
                print(l)
    pex2 = np.hstack((effort[0], effort, effort[-1]))
    pex = effort
    for i in range(matrix.shape[1]):
        aug2 = np.hstack((0, matrix[:, i], 0))
        aug = matrix[:, i]
#        lab = label + " " + str(thresholdlist[i])
        lab = str(thresholdlist[i]) + "x"
        plt.grid(axis='both')
        if SHADED == 0:
            plt.title(label)
            plt.semilogx(pex, aug, "-o", label=lab)
        elif SHADED == 2:
            lab = label + str(thresholdlist[i]) + "x"
            lab = label
            plt.semilogx(pex, aug, "-", label=lab, color=color)
            plt.ylabel("Nonunique fraction of data")
        elif SHADED == 3:
            plt.semilogy(pex, aug, "-", label=lab, color=color)
            plt.ylabel("Number of unique categories ")
            plt.xlabel("Sampling effort")
        elif SHADED == 1:
            plt.subplot(numplots, 1, n + 1)
            plt.semilogx(pex, aug, "-", label=lab, color=color)
            plt.fill(pex2, aug2, "k", alpha=0.2)
            plt.title(label)
            plt.ylabel("Fraction of data")
        else:
            plt.semilogx(pex, aug, "-", label=lab)
            plt.fill(pex2, aug2, "k", alpha=0.2)
            plt.title(label)
            plt.ylabel("Fraction of data")
#            label=str(thresholdlist[i]))
#        plt.fill(pex, aug, "k", alpha=0.2)
    plt.ylim((0, 1))
    plt.xlim((1E4, 1E11))
    if SHADED == 0 or n + 1 == numplots:
        plt.xlabel("Sequencing effort (bp)")
    else:    # suppress drawing of x-axis labels for all but last plot
        frame1 = plt.gca()
        frame1.axes.get_xaxis().set_ticks([])
    plt.tight_layout()
    return()

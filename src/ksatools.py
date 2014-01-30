#!/usr/bin/env python
'''Tool to generate graphs of kmer spectra'''

import sys
import numpy as np

def pad(xvalues, yvalues):
    '''Adds missing integer values to x and corresponding zeros to y.'''
    yout = []
    xout = []
    for i in range(int(1.5*(max(xvalues)))):
        try:
            xout.append(xvalues[xvalues.index(i+1)])
            yout.append(yvalues[xvalues.index(i+1)])
        except ValueError:
            xout.append(i+1)
            yout.append(0)
    return(xout, yout)

def getcolor(index):
    '''returns a string that cycles through more colors than default'''
    colorlist = ["b", "g", "r", "c", "y", "m", "k", "BlueViolet",
            "Coral", "Chartreuse", "DarkGrey", "DeepPink", "LightPink"]
    l = index % len(colorlist)
    return colorlist[l]

def calcmedian(yd, y, num):
    '''interpolates, returning value of yd corresponding to to num on y'''
    try:
        top = np.max(np.nonzero(y > num))
    except ValueError:
        top = None
    try:
        bottom = np.min(np.nonzero(y <= num))
    except ValueError:
        bottom = None
    if top != None and bottom != None:
        cutoff = yd[bottom] + (num - y[bottom]) / (y[top] - y[bottom]) * (yd[top] -yd[bottom])
    elif top != None and bottom == None:
        cutoff = ((num *1.0) / (y[top]) * (yd[top]))
    elif top == None and bottom != None:
        cutoff = yd[bottom]
    else:
        cutoff = 0
    if num <= 1:
        cutoff = np.ceil(cutoff - .001)
    return cutoff

def cleanlabel(label):
    '''Sanitizes graph labels of unintersting file extensions'''
    suffixes = [".histhist", ".fastq", "_info_contigstats.txt",
                ".stats.txt", ".txt", ".csv", ".037.kmerhistogram"]
    for suffix in suffixes:
        if label.find(suffix) > 0:
            label = label[0:(label.find(suffix))]
    return label

def getmgrkmerspectrum(accessionnumber, mgrkey=None):
    '''Retrieve kmer spectrum from MG-RAST'''
    import urllib2, json, time
    assert accessionnumber[0:3] == "mgm", sys.exit(
        "Data error: field %s not in mgm......... accession number format" % accessionnumber)
    some_url = "http://api.metagenomics.anl.gov/api.cgi/metagenome/%s?verbosity=full" % accessionnumber
    if mgrkey != None:
        some_url = some_url + "&auth=%s" % mgrkey
    sys.stderr.write("Sending request for " + some_url + "\n")
    time.sleep(1)
# Ok, exception handling here is a important.  HTTP errors and
# malformed JSON are likely failure modes.
    try:
        opener = urllib2.urlopen(some_url)
    except urllib2.HTTPError, e:
        sys.stderr.write("Error retrieving %s" % some_url)
        sys.stderr.write("Error with HTTP request: %d %s\n%s" % (e.code, e.reason, e.read()))
        return np.atleast_2d(np.array([1, 0]))
    try:
        j = json.loads(opener.read())
    except ValueError:
        sys.stderr.write("Error parsing %s" % some_url)
        j = {}
    try:
        sys.stderr.write("Error with %s\nERROR : %s\n" % (some_url, j["ERROR"]))
        dataarray = None
    except KeyError:
        try:
#        This is the data object containing the 15mer spectrum
            spectrum = j["statistics"]["qc"]["kmer"]["15_mer"]["data"]
            dataarray = np.array(spectrum, dtype="float")
            try:
                dataarray = dataarray[:, 0:2]
            except IndexError:
                dataarray = np.atleast_2d(np.array([1, 0]))
        except KeyError:
            dataarray = np.atleast_2d(np.array([1, 0]))
    return dataarray

def calccumsum(a):
    '''Calcaulates the cumulative-sum vectors from a 2d numpy array
    of [cov, num].  Note depends on upstream sort '''
    cn = a[:, 0]                          #   Coverage
    c1 = a[:, 1]                          #   number of distinct kmers.
    cp = cn * c1  # elementwise multiply     observed kmers by abundance
    yd = np.flipud(np.flipud(c1).cumsum()) # cumul number of distinct kmers (top to bottom)
    yo = np.flipud(np.flipud(cp).cumsum()) # cumul number of observed kmers (top to bottom)
    zd = np.cumsum(c1)                     # cumul number of distinct kmers (bottom to top)
    zo = np.cumsum(cp)                     # cumul number of observed kmers (bottom to top)
    if zo.max() == 0:
        raise ValueError  # There should be data here
    return(cn, c1, yd, yo, zd, zo)

def printstats(a, filename, filehandle=None, n=0):
    '''Prints summary statistics to filename'''
    cn, c1, yd, yo, zd, zo = calccumsum(a)
    T = zo.max()
    j = cn / T
    intermediate = - c1 * j * np.log(j)
    # allow calculation with 0 counts in some rows
    intermediate[np.isnan(intermediate)] = 0
    H = np.exp(sum(intermediate))   # antilog Shannon entropy
    if T == 0:
        H = np.NaN
    H2 = 1 / sum(c1 * j * j)        # antilog Reyni entropy / Simpson index
    w = yo / yo.max()
    wd = yd
    M90 = calcmedian(wd, w, .9)      # 90th percentile by observations
    M50 = calcmedian(wd, w, .5)      # 50th percentile by observations
    M10 = calcmedian(wd, w, .1)      # 10th percentile by observations
    M100 = calcmedian(wd, w, 1.0)    # should be the same as wd.max()
    F100 = calcmedian(w, wd, 100)    # fraction of data in top 100 kmers
    F10K = calcmedian(w, wd, 10000)  # in 10K kmers
    F1M = calcmedian(w, wd, 1000000) # in 1M kmers
    if filehandle == None:
        consensusfh = open(filename, "w")
    else:
        consensusfh = filehandle
    if filehandle == None or n == 0:
        consensusfh.write("#filename\tM10\tM50\tM90\tM100\tF100\tF10K" +
             "\tF1M\tH\tH2\n")
    consensusfh.write("%s\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%.1f\t%.1f\n" %
                      (filename, M10, M50, M90, M100, F100, F10K, F1M, H, H2))
    if filehandle == None:
        consensusfh.close()

def loadfile(filename):
    '''Loads file, returns two-column ndarray or None on
    failure.  Uses filename to guess format.'''
    try:
        # parse velvet contig stats format
        if filename.find("stats.txt") >= 0:
            matrix = np.loadtxt(filename, usecols=(5, 1), skiprows=1)
        else: # default bare-bones spectrum format
            matrix = np.loadtxt(filename)
        # return None if the file is empty
        if matrix.shape[0] == 0:
            return None
        else:
            return np.atleast_2d(matrix)
    except IOError:
        sys.stderr.write("ERROR: Can't find file %s\n" % filename)
        return None

def makegraphs(spectrum, filename, option=6, label=None, n=0, dump=False):
    '''Draw graphs, one at a time, and add them to the current plot.
    spectrum contains the data; filename is the file stem for saving
    option determines the type of graph; label labels each trace;
    n counts the (successful) traces.  Returns n.'''
    import matplotlib.pyplot as plt
    # note, calccumsum will raise an exception here if data is invalid
    (cn, c1, yd, yo, zd, zo) = calccumsum(spectrum)
    if label == None:
        tracelabel = cleanlabel(filename)
    else:
        tracelabel = cleanlabel(label)
    assert spectrum.dtype == "float"
    # sorted by abundance/coverage
    b = np.flipud(spectrum[np.argsort(spectrum[:, 0]), :])
    # sorted by size (for contigs)
    c = np.flipud(spectrum[np.argsort(spectrum[:, 1]), :])
    # sorted by abundance-size product (explained)
    d = np.flipud(spectrum[np.argsort(spectrum[:, 1] * spectrum[:, 0]), :])
    (b_cn, b_c1, b_yd, b_yo, b_zd, b_zo) = calccumsum(b) # abundance
    (c_cn, c_c1, c_yd, c_yo, c_zd, c_zo) = calccumsum(c) # size
    (d_cn, d_c1, d_yd, d_yo, d_zd, d_zo) = calccumsum(d) # explained
    No = b_zo.max()
    Nd = b_zd.max()
    x = np.arange(len(b[:, 0]))                  # rank
    color = getcolor(n)
    if option == 0:
        pA = plt.loglog(b_cn, b_c1, "-", color=color, label=tracelabel)
        pA = plt.loglog(b_cn, b_c1, ".", color=color)
        plt.xlabel("kmer abundance")
        plt.ylabel("number of kmers")
        legendloc = "upper right"
        plt.grid(1)
    if option == 0 or option == -1:
        if dump:
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
        legendloc = "upper right"
        plt.grid(1)
        if dump:
            c = np.hstack((cn.reshape((len(cn), 1)),
                ((cn * c1).reshape((len(cn), 1)))))
            sys.stderr.write("saving output table in %s.1.plot.csv\n" %
                filename)
            np.savetxt("%s.1.plot.csv" % filename, c, fmt=['%d', '%d'],
                delimiter="\t")
    elif option == 2:
        pA = plt.loglog(b_zo, b_cn, color=color, label=tracelabel)
        plt.xlabel("fraction of data")  # cumulative kmers observed
        plt.ylabel("kmer abundance")
        legendloc = "lower left"
        plt.grid(1)
    elif option == 3:
        pA = plt.semilogy(b_zo / No, b_cn, color=color, label=tracelabel)
        pA = plt.semilogy(b_zo / No, b_cn, '.', color=color)
        plt.xlabel("fraction of data")  # formerly "fraction of observed kmers"
        plt.ylabel("kmer abundance ")
        plt.grid(1)
        legendloc = "lower left"
    elif option == 4: # Fraction of distinct kmers vs abundance  NOT RECOMMENDED
        pA = plt.semilogy(b_zd / Nd, b_cn, color=color, label=tracelabel)
        pA = plt.semilogy(b_zd / Nd, b_cn, '.', color=color)
        plt.xlabel("fraction of distinct kmers")
        plt.ylabel("kmer abundance")
        legendloc = "upper right"
        plt.grid(1)
    elif option == 5:
        pA = plt.semilogx(yd, zo / No, '-', color=color)
        pA = plt.semilogx(yd, zo / No, '.', color=color, label=tracelabel)
        plt.xlabel("kmer rank (bp)")
        plt.ylabel("fraction of data")
        plt.xlim((1, 10**10))
        plt.ylim(0, 1)
        plt.grid(1)
        legendloc = "lower left"
    elif option == 6:
        pA = plt.loglog(b_zd, b_cn, '-', color=color, label=tracelabel)
        pA = plt.loglog(b_zd, b_cn, '.', color=color)
        plt.xlabel("kmer rank")
        plt.ylabel("kmer abundance")
        plt.xlim((1, 10**10))
        plt.ylim(1, 10**7)
        plt.grid(1)
        legendloc = "lower left"
        if dump:
            c = np.hstack((yd.reshape((len(yd), 1)), cn.reshape((len(cn), 1))))
            sys.stderr.write("saving output table in %s.6.plot.csv\n" % filename)
            np.savetxt("%s.6.plot.csv" % filename, c, fmt=['%d', '%d'], delimiter="\t")

    elif option == 7:
        pA = plt.plot(x, c_zd, '.-', color=color, label=tracelabel)
        plt.xlabel("contig size rank")
        plt.ylabel("cuml contig size")
        plt.grid(1)
        legendloc = "upper right"
    elif option == 8:
        pA = plt.plot(x, c_zo / No, '.-', color=color, label=tracelabel)
        plt.xlabel("contig size rank ")
        plt.ylabel("frac data explained ")
        plt.grid(1)
        legendloc = "upper right"
    elif option == 9:
        pA = plt.plot(x, d_zo, '-', color=color, label=tracelabel)
        plt.xlabel("contig explain rank ")
        plt.ylabel("data explained ")
        plt.grid(1)
        legendloc = "upper right"
    elif option == 10:
        pA = plt.plot(x, b_yo / No, '.-', color=color, label=tracelabel)
        plt.xlabel("contig cov rank ")
        plt.ylabel("frac data explained ")
        plt.grid(1)
        legendloc = "upper right"
    elif option == 11:
        pA = plt.plot(x, b_yo, '.-', color=color, label=tracelabel)
        plt.xlabel("contig cov rank ")
        plt.ylabel("data explained (bogo bp) ")
        plt.grid(1)
        legendloc = "upper right"
    elif option == 12:
        pA = plt.plot(c_zd, c_zo, '.-', color=color, label=tracelabel)
        plt.xlabel("cumulative contig size")
        plt.ylabel("data explained (bogo bp) ")
        plt.grid(1)
        legendloc = "upper right"
    elif option == 13:
        pA = plt.plot(x, b_cn, '.-', color=color, label=tracelabel)
        plt.xlabel("contig cov rank ")
        plt.ylabel("kmer abundance")
        plt.grid(1)
        legendloc = "upper right"
    elif option == 14:
        pA = plt.plot(x, d_cn * d_c1, '.-', color=color, label=tracelabel)
        plt.xlabel("contig expl rank ")
        plt.ylabel("data explained (bogo bp) ")
        plt.grid(1)
        legendloc = "upper right"
    elif option == 15:
        pA = plt.plot(x, c_c1, '.-', color=color, label=tracelabel)
        plt.xlabel("contig size rank")
        plt.ylabel("contig size")
        plt.grid(1)
        legendloc = "upper right"
    elif option == 16:
        pA = plt.plot(x, c_yd, '.-', color=color, label=tracelabel)
        plt.xlabel("contig expl rank ")
        plt.ylabel("data explained (bogo bp) ")
        plt.grid(1)
        legendloc = "upper right"
    elif option == 17:
        pA = plt.plot(x, d_cn * d_c1, '.-', color=color, label=tracelabel)
        plt.xlabel("contig expl rank ")
        plt.ylabel("data explained (bogo bp) ")
        plt.grid(1)
        legendloc = "upper right"
    if option >= 0:
        plt.legend(loc=legendloc)

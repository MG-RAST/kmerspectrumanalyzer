#!/usr/bin/env python
'''Tool to generate graphs of kmer spectra'''

import sys
import numpy as np

def renyispectrum(x, spectrum):
    n = spectrum[:, 0]
    y = spectrum[:, 1]
    N = np.sum(n*y)
    p = n / N
    R = np.zeros(x.shape)
    for i, l in enumerate(x):
        G = np.log10(np.sum(y * np.power(p, l)))/(1-l)
        R[i] = G
    return R

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
    return cutoff

def cleanlabel(label):
    '''Sanitizes graph labels of unintersting file extensions'''
    suffixes = [".histhist", ".fastq", "_info_contigstats.txt",
                ".stats.txt", ".txt", ".csv", ".037.kmerhistogram"]
    for suffix in suffixes:
        if label.find(suffix) > 0:
            label = label[0:(label.find(suffix))]
    return label

def drawboxes(breaks, axis, boxcolor=1):
    '''Draws boxes on the current plot.'''
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    import matplotlib.pyplot as plt
    ax = plt.gca()
    x1, x2 = plt.xlim()
    y1, y2 = plt.ylim()
    patches = []
    if axis == 0:
        for i in range(len(breaks)-1):
            y1, y2 = (breaks[i+1], breaks[i])
            patches.append(Polygon([[x1, y2], [x1, y1], [x2, y1], [x2, y2]], True))
    else:
        for i in range(len(breaks)-1):
            x1, x2 = (breaks[i+1], breaks[i])
            patches.append(Polygon([[x1, y2], [x1, y1], [x2, y1], [x2, y2]], True))
    if boxcolor == 1:
        p = PatchCollection(patches, cmap=plt.cm.jet, alpha=0.4)
    else:
        p = PatchCollection(patches, cmap=plt.cm.Greys, alpha=0.2)
    p.set_array(np.array([0, 0.2, 0.4, 0.5, 0.7, 0.9, 1]))
    ax.add_collection(p)

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
    w = np.array(yo, dtype="float") / yo.max()
    wd = yd
    M90 = calcmedian(wd, w, .9)      # 90th percentile by observations
    M50 = calcmedian(wd, w, .5)      # 50th percentile by observations
    M10 = calcmedian(wd, w, .1)      # 10th percentile by observations
    M100 = calcmedian(wd, w, 1.0)    # should be the same as wd.max()
    F100 = calcmedian(w, wd, 100)    # fraction of data in top 100 kmers
    F10K = calcmedian(w, wd, 10000)  # in 10K kmers
    F1M = calcmedian(w, wd, 1000000) # in 1M kmers
    C50 = calcmedian(cn, w, .5)      # 50th percentile by observations
    AVONE = np.sum(cn * c1) / T
    assert (AVONE - 1.) < 1E-7
    AVCOV = np.sum(cn * c1 * cn) / T
    AVGEN = np.sum(cn * c1 * c1) / T
    if filehandle == None:
        consensusfh = open(filename, "w")
    else:
        consensusfh = filehandle
    if filehandle == None or n == 0:
        consensusfh.write("#filename\tM10\tM50\tM90\tM100\tF100\tF10K" +
             "\tF1M\tH\tH2\tAVC\tAVG\tC50\n")
    consensusfh.write("%s\t%.1f\t%.1f\t%.1f\t%d\t%f\t%f\t%f\t%.1f\t%.1f\t%.1f\t%.1f\t%.2f\n" %
                      (filename, M10, M50, M90, M100, F100, F10K, F1M, H, H2, AVCOV, AVGEN, C50))
    if filehandle == None:
        consensusfh.close()

def loadfile(filename):
    '''Loads file, returns two-column ndarray or None on
    failure.  Uses filename to guess format.'''
    try:
        # parse velvet contig stats format
        if filename.find("stats.txt") >= 0:
            matrix = np.loadtxt(filename, usecols=(5, 1), skiprows=1)
        if filename.find(".npo") >= 0:
            matrix = np.loadtxt(filename, usecols=(0, 1), skiprows=6)
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

def plotstratify(spectrum, bands=None):
    bands, frac, size = stratify(spectrum, bands)
    for i in range(len(bands)):
        if i != len(bands)-1:
            print "%.04f"%(frac[i] - frac[i+1]), "\t% 13d"% (
                size[i] - size[i+1]), "\t", str(bands[i])+"-"+str(bands[i+1])

def stratify(spectrum, bands=None):
    '''Breaks spectrum up into defined abundance-buckets,
    reporting data fraction and number of kmers=basepairs
    contained in each bucket.'''
    if bands == None:
        bands = [1, 10, 100, 1000, 10000, 100000]
    n = spectrum[:, 0]
    y = spectrum[:, 1]
    p = n * y
    T = sum(p)
    frac = []
    size = []
    for b in bands:
        frac.append(np.sum(p[n >= b]) / T)
        size.append(np.sum(y[n >= b]) + 1)
    # +1 adds an artificial count in the largest bin
    frac.append(0)
    size.append(0)
    bands.append(bands[-1] * 10)
    return bands, frac, size

def makegraphs(spectrum, filename, option=6, label=None, n=0, dump=False, opts=None):
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
    outfile = "%s.%d.plot.csv" % (filename, option)
    if option == 0:
        pA = plt.loglog(b_cn, b_c1, "-", color=color, label=tracelabel)
        pA = plt.loglog(b_cn, b_c1, ".", color=color)
        xlabel, ylabel = ("kmer abundance", "number of kmers")
        legendloc = "upper right"
    if option == 0 or option == -1:
        if dump:
            c = np.hstack((cn.reshape((len(cn), 1)),
                (c1.reshape((len(cn), 1)))))
            sys.stderr.write("saving output table in %s\n" % outfile)
            np.savetxt(outfile, c, fmt=['%d', '%d'], delimiter="\t")
    elif option == 1:
        pA = plt.loglog(cn, cn * c1, "-", color=color, label=tracelabel)
        pA = plt.loglog(cn, cn * c1, '.', color=color)
        xlabel, ylabel = ("kmer abundance", "kmers observed")
        legendloc = "upper right"
        if dump:
            c = np.hstack((cn.reshape((len(cn), 1)),
                ((cn * c1).reshape((len(cn), 1)))))
            sys.stderr.write("saving output table in %s\n" % outfile)
            np.savetxt(outfile, c, fmt=['%d', '%d'],
                delimiter="\t")
    elif option == 2:
        pA = plt.loglog(b_zo, b_cn, color=color, label=tracelabel)
        xlabel, ylabel = ("basepairs observed", "kmer abundance")
        legendloc = "lower left"
    elif option == 3:
        pA = plt.semilogy(b_zo / No, b_cn, color=color, label=tracelabel)
        pA = plt.semilogy(b_zo / No, b_cn, '.', color=color)
        xlabel, ylabel = ("fraction of observed data", "kmer abundance")
        legendloc = "lower left"
    elif option == 4: # Fraction of distinct kmers vs abundance  NOT RECOMMENDED
        pA = plt.semilogy(b_zd / Nd, b_cn, color=color, label=tracelabel)
        pA = plt.semilogy(b_zd / Nd, b_cn, '.', color=color)
        xlabel, ylabel = ("fraction of distinct kmers", "kmer abundance")
        legendloc = "upper right"
    elif option == 5 or option == 25 or option == 24:
        pA = plt.semilogx(yd, zo / No, '-', color=color)
        pA = plt.semilogx(yd, zo / No, '.', color=color, label=tracelabel)
        xlabel, ylabel = ("kmer rank (bp)", "fraction of observed data")
        plt.xlim((1, 10**10))
        plt.ylim(0, 1)
        legendloc = "lower left"
    elif option == 6 or option == 26:
        pA = plt.loglog(b_zd, b_cn, '-', color=color, label=tracelabel)
        pA = plt.loglog(b_zd, b_cn, '.', color=color)
        xlabel, ylabel = ("kmer rank (bp)", "kmer abundance")
        plt.xlim((1, 10**10))
        if max(b_cn) < 10**8:
            plt.ylim(1, 10**7)
        legendloc = "lower left"
        if dump:
            c = np.hstack((yd.reshape((len(yd), 1)), cn.reshape((len(cn), 1))))
            sys.stderr.write("saving output table in %s\n" % outfile)
            np.savetxt(outfile, c, fmt=['%d', '%d'], delimiter="\t")
    elif option == 7:
        pA = plt.plot(x, c_zd, '.-', color=color, label=tracelabel)
        xlabel, ylabel = ("contig size rank", "cuml contig size")
        legendloc = "upper right"
    elif option == 8:
        pA = plt.plot(x, c_zo / No, '.-', color=color, label=tracelabel)
        xlabel, ylabel = ("contig size rank", "frac data explained")
        legendloc = "upper right"
    elif option == 9:
        pA = plt.plot(x, d_zo, '-', color=color, label=tracelabel)
        xlabel, ylabel = ("contig explain rank", "data explained")
        legendloc = "upper right"
    elif option == 10:
        pA = plt.plot(x, b_yo / No, '.-', color=color, label=tracelabel)
        xlabel, ylabel = ("contig cov rank", "frac data explained")
        legendloc = "upper right"
    elif option == 11:
        pA = plt.plot(x, b_yo, '.-', color=color, label=tracelabel)
        xlabel, ylabel = ("contig cov rank", "data explained (bogo bp)")
        legendloc = "upper right"
    elif option == 12:
        pA = plt.plot(c_zd, c_zo, '.-', color=color, label=tracelabel)
        xlabel, ylabel = ("cumulative contig size", "data explained (bogo bp)")
        legendloc = "upper right"
    elif option == 13:
        pA = plt.plot(x, b_cn, '.-', color=color, label=tracelabel)
        xlabel, ylabel = ("contig cov rank", "kmer abundance")
        legendloc = "upper right"
    elif option == 14:
        pA = plt.plot(x, d_cn * d_c1, '.-', color=color, label=tracelabel)
        xlabel, ylabel = ("contig expl rank", "data explained (bogo bp)")
        legendloc = "upper right"
    elif option == 15:
        pA = plt.plot(x, c_c1, '.-', color=color, label=tracelabel)
        xlabel, ylabel = ("contig size rank", "contig size")
        legendloc = "upper right"
    elif option == 16:
        pA = plt.plot(x, c_yd, '.-', color=color, label=tracelabel)
        xlabel, ylabel = ("contig expl rank", "data explained (bogo bp)")
        legendloc = "upper right"
    elif option == 17:
        pA = plt.plot(x, d_cn * d_c1, '.-', color=color, label=tracelabel)
        xlabel, ylabel = ("contig expl rank", "data explained (bogo bp)")
        legendloc = "upper right"
    elif option == 18:  # semilog version of 1
        pA = plt.semilogx(b_cn, b_c1, '.-', color=color, label=tracelabel)
        xlabel, ylabel = ("number of reads sampled", "fraction nonunique")
        legendloc = "lower right"
    elif option == 30:
        lam = np.arange(.01, 10, .01)
        entropyspectrum = np.power(10,renyispectrum(lam, spectrum))
        pA = plt.semilogy(lam, entropyspectrum, '.-', color=color, label=tracelabel)
        xlabel, ylabel = ("lambda", "Renyientropy")
        legendloc = "upper right"

    if option == -2 or option == 26 or option == 25:
        if dump:
            plotstratify(spectrum)
        else:
            plotstratify(spectrum)
    # For these two graphs, draw rainbow bands
    if option == 26:
        bands, frac, size = stratify(spectrum)
        drawboxes(bands, 0)
    elif option == 25:
        bands, frac, size = stratify(spectrum)
        fracboundaries = 1 - np.array(frac)
        sizeboundaries = size
        drawboxes(sizeboundaries, 1)
        drawboxes(fracboundaries, 0, boxcolor=0)
    elif option == 24:
        bands, frac, size = stratify(spectrum)
        fracboundaries = 1 - np.array(frac)
        sizeboundaries = size
        drawboxes(fracboundaries, 0)

    # Draw graphs if option >= 0
    if option >= 0:
        if not opts.suppress:
            plt.legend(loc=legendloc)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if hasattr(opts, "title") and not opts.title == None and n == 0:
            plt.title(opts.title)
        plt.grid(1)

#!/usr/bin/env python
from __future__ import division

'''Tool to generate graphs of kmer spectra'''

COLORLIST = [
    "b", "g", "r", "c", "y",
    "m", "k", "BlueViolet", "Coral", "Chartreuse",
    "DarkGrey", "DeepPink", "LightPink"]

import sys
import numpy as np

def renyispectrum(x, spectrum):
    '''Given a two-dimensional spectrum spectrum and a one-dimensional
    vector of Renyi spectrum exponents, return a vector of corresponding
    log10 Renyi entropies.
    '''
    n = spectrum[:, 0]
    y = spectrum[:, 1]
    N = np.sum(n*y)
    p = n / N
    R = np.zeros(x.shape)
    for i, l in enumerate(x):
        G = np.log10(np.sum(y * np.power(p, l)))/(1-l)
        R[i] = G
    return R

def pad(xvalues, yvalues, fill=True):
    '''Adds missing integer values to x and corresponding zeros to y.'''
    yout = []
    xout = []
    xv = np.arange(1, int(1.5*max(xvalues)))
    yv = np.zeros(shape=int(1.5*max(xvalues)-1), dtype=int)

    if fill == True:
        for i, x in enumerate(xvalues):
            yv[np.where(xv == x)[0][0]] = yvalues[i]
        xout = list(xv)
        yout = list(yv)
    else:
        for i, x in enumerate(xvalues):
            if (x-1 not in xvalues) and (x-1) not in xout:
                xout.append(x-1)
                yout.append(0)
            if  yvalues[xvalues.index(x)] > 0:
                xout.append(xvalues[xvalues.index(x)])
                yout.append(yvalues[xvalues.index(x)])
            if (x+1 not in xvalues) and (x+1) not in xout:
                xout.append(x+1)
                yout.append(0)
    return(xout, yout)

def smoothspectrum(data):
    bins = np.hstack((np.arange(1, 200), np.exp(np.arange(0, 4, .01) * np.log(10)) * 200))
    x = data[:, 0]
    y = data[:, 1]
    yo = np.zeros((len(bins)))
    for i in range(0, len(bins)-1):
        yo[i] = y[np.where((x >= bins[i])* (x < bins[i+1]))].sum()
    xo = bins[:-1]
    yo = yo[:-1] / np.diff(bins)
    return xo, yo

def getcolor(index, colorlist):
    if colorlist == []:
        colorlist = COLORLIST
    l = index % len(colorlist)
    return colorlist[l]

def calcmedian(yd, y, num):
    '''wrapper for np.interp; interpolates to return the value of yd corresponding to num on y
    sorts the data and prepends a 0,0 to get smooth behavior for entire range 0,1.'''
    ya = np.argsort(y)
    y2 = np.hstack(([0], y[ya]))
    y2d = np.hstack(([0], yd[ya]))
    r = np.interp(num, y2, y2d)
    return r

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
    except urllib2.HTTPError as e:
        sys.stderr.write("Error retrieving %s" % some_url)
        sys.stderr.write("Error with HTTP request: %d %s\n%s" % (e.code, e.reason, e.read()))
        return np.atleast_2d(np.array([1, 0]))
    try:
        j = json.loads(opener.read())
    except ValueError:
        sys.stderr.write("Error parsing %s" % some_url)
        j = {}
    try:
        sys.stderr.write("Error with %s\nERROR: %s\n" % (some_url, j["ERROR"]))
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
    if filehandle is None:
        consensusfh = open(filename, "w")
    else:
        consensusfh = filehandle
    if filehandle is None or n == 0:
        consensusfh.write(
            "#filename\tM10\tM50\tM90\tM100\tF100\tF10K" +
            "\tF1M\tH\tH2\tAVC\tAVG\tC50\n")
    consensusfh.write(
        "%s\t%.1f\t%.1f\t%.1f\t%d\t%f\t%f\t%f\t%.1f\t%.1f\t%.1f\t%.1f\t%.2f\n" %
        (filename, M10, M50, M90, M100, F100, F10K, F1M, H, H2, AVCOV, AVGEN, C50))
    if filehandle is None:
        consensusfh.close()

def getlength(filename):
    '''parse nonpareil output npo for length parameter'''
    for line in open(filename):
        if line[0:5] == "# @L:":
            return float(line[6:])

def loadfile(filename):
    '''Loads file, returns two-column ndarray or None on
    failure.  Uses filename to guess format.'''
    try:
        # parse velvet contig stats format
        if filename.find("stats.txt") >= 0:
            matrix = np.loadtxt(filename, usecols=(5, 1), skiprows=1)
        elif filename.find(".npo") == len(filename)-4:
            matrix = np.loadtxt(filename, usecols=(0, 1), skiprows=6)
            L = getlength(filename)
            matrix[:, 0] = matrix[:, 0] * L
        else: # default bare-bones spectrum format
            try:
                matrix = np.loadtxt(filename, comments="#")
            except ValueError:
                matrix = np.loadtxt(filename, skiprows=1, delimiter=",", usecols=(0, 1))
        # return None if the file is empty
        matrix[np.isinf(matrix)] = 0
        matrix = np.atleast_2d(matrix)
        if matrix.shape[0] == 0:
            return []
        else:
            return matrix
    except IOError:
        sys.stderr.write("ERROR: Can't find file %s\n" % filename)
        return []

def printstratify(spectrum, bands=None, flat=False, label=""):
    '''prints table of sizes and data fracitons for separate
    bands, inclusive of lower bound, excluding upper band
    boundary'''
    bands, frac, size = stratify(spectrum, bands)
    if flat is False:
        for i in range(len(bands)):
            if i != len(bands)-1:
                print("{:.04f}\t{: 13.0f}\t{:.0f}-{:.0f}\n".format((frac[i] - frac[i+1]), (
                    size[i] - size[i+1]), bands[i], bands[i+1]))
    else:
        print("#name\t"+"\t".join(map(str, list(bands[:-1]+bands[1:]))))
        print(label+"\t"+"\t".join(map(str, list(size)[:-1] +
                                       list(frac)[1:])))

def stratify(spectrum, bands=None):
    '''Breaks spectrum up into defined abundance-buckets,
    reporting cumulative data fraction and number of
    kmers=basepairs with kmer-depths greater than or equal
    to each bucket boundary.  Returns list of bands, list
    of (cumulative) fractions, and list of (cumulative)
    basepairs for each band.
    '''
    if bands is None:
        bands = [1, 10, 100, 1000, 10000, 100000]
    n = spectrum[:, 0]
    y = spectrum[:, 1]
    p = n * y
    T = sum(p)
    frac = []
    size = []
    for b in bands:
        frac.append(np.sum(p[n >= b]) / T)
        size.append(np.sum(y[n >= b]))
    return bands, frac, size

def makegraphs(spectrum, filename, option=6, label=None, n=0,
               dump=False, opts=None, colorlist=COLORLIST, 
               stylelist=None):
    '''Draw graphs, one at a time, and add them to the current plot.
    spectrum contains the data; filename is the file stem for saving
    option determines the type of graph; label labels each trace;
    n counts the (successful) traces.  Returns n.
    '''
    import matplotlib.pyplot as plt
    # note, calccumsum will raise an exception here if data is invalid
    if label is None:
        tracelabel = cleanlabel(filename)
    else:
        tracelabel = cleanlabel(label)
    if option == 0 or option == 1 or option == 19 or option == 20:
        xclean, yclean = smoothspectrum(spectrum)
        spectrum = np.vstack([xclean, yclean]).T
    if option == 21:
        xclean, yclean = pad(list(spectrum[:, 0]), list(spectrum[:, 1]), fill=False)
        spectrum = np.vstack([xclean, yclean]).T
    assert spectrum.dtype == "float"
    (cn, c1, yd, yo, zd, zo) = calccumsum(spectrum)
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
    color = getcolor(n, colorlist)
    outfile = "%s.%d.plot.csv" % (filename, option)
    style = ".-"
    drawstyle = None
    if option == 0 or True:  # fallback visualization
        plot1, p, q = (plt.loglog, b_cn, b_c1)
        xlabel, ylabel = ("kmer abundance", "number of kmers")
        legendloc = "upper right"
    if option == 1:
        plot1, p, q = (plt.loglog, cn, cn * c1)
        xlabel, ylabel = ("kmer abundance", "kmers observed")
        legendloc = "upper right"
    elif option == 2:
        plot1, p, q = (plt.loglog, b_zo, b_cn)
        xlabel, ylabel = ("basepairs observed", "kmer abundance")
        legendloc = "lower left"
    elif option == 3:
        plot1, p, q = (plt.semilogy, b_zo / No, b_cn)
        xlabel, ylabel = ("fraction of observed data", "kmer abundance")
        legendloc = "lower left"
    elif option == 4: # Fraction of distinct kmers vs abundance  NOT RECOMMENDED
        plot1, p, q = (plt.semilogy, b_zd / Nd, b_cn)
        xlabel, ylabel = ("fraction of distinct kmers", "kmer abundance")
        legendloc = "upper right"
    elif option == 5 or option == 25 or option == 24:
        plot1, p, q = (plt.semilogx, yd, 1-zo/No)
        xlabel, ylabel = ("kmer rank (bp)", "cumulative fraction of observed data")
        plt.xlim((1, 10**11))
        plt.ylim(0, 1)
        legendloc = "lower left"
    elif option == 6 or option == 26:
        plot1, p, q = (plt.loglog, b_zd, b_cn)
        xlabel, ylabel = ("kmer rank (bp)", "kmer abundance")
        plt.xlim((1, 10**11))
        if max(b_cn) < 10**8:
            plt.ylim(1, 10**7)
        legendloc = "lower left"
    elif option == 7:
        plot1, p, q = (plt.plot, x, c_zd)
        xlabel, ylabel = ("contig size rank", "cuml contig size")
        legendloc = "upper right"
    elif option == 8:
        plot1, p, q = (plt.plot, x, c_zo / No)
        xlabel, ylabel = ("contig size rank", "frac data explained")
        legendloc = "upper right"
    elif option == 9:
        plot1, p, q = (plt.plot, x, d_zo)
        xlabel, ylabel = ("contig explain rank", "data explained")
        legendloc = "upper right"
    elif option == 10:
        plot1, p, q = (plt.plot, x, b_yo / No)
        xlabel, ylabel = ("contig cov rank", "frac data explained")
        legendloc = "upper right"
    elif option == 11:
        plot1, p, q = (plt.plot, x, b_yo)
        xlabel, ylabel = ("contig cov rank", "data explained (bogo bp)")
        legendloc = "upper right"
    elif option == 12:
        plot1, p, q = (plt.plot, c_zd, c_zo)
        xlabel, ylabel = ("cumulative contig size", "data explained (bogo bp)")
        legendloc = "upper right"
    elif option == 13:
        plot1, p, q = (plt.plot, x, b_cn)
        xlabel, ylabel = ("contig cov rank", "kmer abundance")
        legendloc = "upper right"
    elif option == 14:
        plot1, p, q = (plt.plot, x, d_cn * d_c1)
        xlabel, ylabel = ("contig expl rank", "data explained (bogo bp)")
        legendloc = "upper right"
    elif option == 15:
        plot1, p, q = (plt.plot, x, c_c1)
        xlabel, ylabel = ("contig size rank", "contig size")
        legendloc = "upper right"
    elif option == 16:
        plot1, p, q = (plt.plot, x, c_yd)
        xlabel, ylabel = ("contig expl rank", "data explained (bogo bp)")
        legendloc = "upper right"
    elif option == 17:
        plot1, p, q = (plt.plot, x, d_cn * d_c1)
        xlabel, ylabel = ("contig expl rank", "data explained (bogo bp)")
        legendloc = "upper right"
    elif option == 18:  # semilog version of 1
        plot1, p, q = (plt.semilogx, b_cn, b_c1)
        xlabel, ylabel = ("number of bases sampled", "nonunique fraction")
        plt.ylim(0, 1)
        legendloc = "lower left"
    elif option == 19:  # stairstep version of 1
        plot1, p, q = (plt.loglog, cn, cn * c1)
        style, drawstyle = ("-", "steps")
        xlabel, ylabel = ("kmer abundance", "kmers observed")
        legendloc = "upper right"
    elif option == 20:  # stairstep version of 18
        plot1, p, q = (plt.semilogx, cn, cn * c1 * cn)
        style, drawstyle = ("-", "steps-mid")
        xlabel, ylabel = ("kmer abundance", "data quantity")
        legendloc = "upper right"
    elif option == 21:  # stairstep, straight axes version of 1
        plot1, p, q = (plt.plot, cn, cn * c1)
        style, drawstyle = ("-", "steps-mid")
        xlabel, ylabel = ("kmer abundance", "kmers observed")
        legendloc = "upper right"
    elif option == 27:
        fminusparameter = np.exp(np.diff(np.log(c1)) - np.log(cn)[1:])
        plot1, p, q = (plt.semilogy, cn[:-1], fminusparameter)
        xlabel, ylabel = ("kmer abundnace", "fminus paraemter")
        legendloc = "upper right"
    elif option == 28:
        plot1, p, q = (plt.plot, b_zo / No, b_cn)
        xlabel, ylabel = ("fraction of observed data", "kmer abundance")
        legendloc = "lower left"
    elif option == 30:
        lam = np.arange(.01, 10, .01)
        entropyspectrum = np.power(10, renyispectrum(lam, spectrum))
        plot1, p, q = (plt.semilogy, lam, entropyspectrum)
        xlabel, ylabel = ("lambda", "Renyi entropy")
        legendloc = "upper right"
    if dump:
        c = np.hstack((p.reshape((len(p), 1)), q.reshape((len(p), 1))))
        sys.stderr.write("saving output table in %s\n" % outfile)
        ptype, qtype = ("%f", "%f")
        if min(p) == 1:
            ptype = "%d"
        if min(q) == 1:
            qtype = "%d"
        np.savetxt(outfile, c, fmt=[ptype, qtype], delimiter="\t",
                   header=xlabel+"\t"+ylabel)

    if option == -2 or option == 26 or option == 25:
        if dump:
            printstratify(spectrum)
        else:
            printstratify(spectrum)
    if option == -3:
        printstratify(spectrum, flat=True, label=tracelabel)
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
    if hasattr(opts, "xlabel") and opts.xlabel is not None:
        xlabel = opts.xlabel
    if hasattr(opts, "ylabel") and opts.ylabel is not None:
        ylabel = opts.ylabel
    if hasattr(opts, "suppress"):
        suppress = opts.suppress
    else:
        suppress = False
    # Draw graphs if option >= 0
    if stylelist is not None and stylelist != []:
        style = stylelist[n]
    if option >= 0:
        plot1(p, q, style, color=color, label=tracelabel, drawstyle=drawstyle)
        if not suppress:
            plt.legend(loc=legendloc)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if hasattr(opts, "title") and not opts.title is None and n == 0:
            plt.title(opts.title)
        plt.grid(1)

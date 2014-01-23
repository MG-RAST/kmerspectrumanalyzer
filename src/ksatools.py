#!/usr/bin/env python
'''Tool to generate graphs of kmer spectra'''

import sys
import numpy as np

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
    except:
        top = None
    try:
        bottom = np.min(np.nonzero(y <= num))
    except:
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

def getmgrkmerspectrum(accessionnumber, MGRKEY=None):
    '''Retrieve kmer spectrum from MG-RAST'''
    import urllib2, json, time
    assert accessionnumber[0:3] == "mgm", sys.exit("Data error: field %s not in mgm......... accession number format" % accessionnumber)
    some_url = "http://api.metagenomics.anl.gov/api.cgi/metagenome/%s?verbosity=full" % accessionnumber
    if MGRKEY != None:
        some_url = some_url + "&auth=%s" % MGRKEY
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

def sortbycp(data):
    CP = np.concatenate((data, np.atleast_2d(data[:, 0] * data[:, 1]).T), axis=1)
    S = []
    for c in np.argsort(CP[:, 2]):
        S.append(data[c, :])
    A = (np.flipud(np.array(S)))
    return A

def calccumsum(a):
    '''Calcaulates the cumulative-sum vectors from a 2d numpy array
    of [cov, num].  Note depends on upstream sort '''
    cn = a[:, 0]                          #   Coverage
    c1 = a[:, 1]                          #   number of distinct kmers.
    cp = cn * c1  # elementwise multiply     observed kmers by abundance
    yd = np.flipud(np.flipud(c1).cumsum()) # cumulative number of distinct kmers (top to bottom)
    yo = np.flipud(np.flipud(cp).cumsum()) # cumulative number of observed kmers (top to bottom)
    zd = np.cumsum(c1)                     # cumulative number of distinct kmers (bottom to top)
    zo = np.cumsum(cp)                     # cumulative number of observed kmers (bottom to top)
    if zo.max() == 0:
        raise Exception  # There should be data here
    y = zo / zo.max()
    return(cn, c1, yd, yo, zd, zo, y)

def printstats(a, filename, filehandle=None, n=0):
    '''Prints summary statistics to filename'''
    cn, c1, yd, yo, zd, zo, y = calccumsum(a)
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
        consensusfh.write("#filename\tM10\tM50\tM90\tM100\tF100\tF10K\tF1M\tH\tH2\n")
    consensusfh.write("%s\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%.1f\t%.1f\n" %
                      (filename, M10, M50, M90, M100, F100, F10K, F1M, H, H2))
    if filehandle == None:
        consensusfh.close()

def loadfile(filename):
    '''Loads file, returns two-column ndarray or None'''
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


#!/usr/bin/env python
'''This script plots 2d histograms'''

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import argparse


def decimate(a, d):
    '''subroutine to reduce the dimensionality of a rectangular (2d histogram) matrix
     parameters are the original matrix a and a number or tuple with the decimation factor'''
    try:
        d2 = d[1]
        d1 = d[0]
    except:
        d2 = d
        d1 = d
    print("D1 %d D2 %d" % (int(d1), int(d2)))
    r1 = a.shape[0]
    r2 = a.shape[1]
    r1s = int(r1 / d1)
    r2s = int(r2 / d2)
    b = []
    for i in range(0, r1s):
        b.append([])
        for j in range(0, r2s):
            s = 0
            for l in range(0, d1):
                for m in range(0, d2):
                    s = s + a[i * d1 + l, j * d2 + m]
            b[i].append(s)
    print("decimate return:", np.array(b).shape)
    return np.array(b)


def decimate_withcovest(a, covest):
    '''attempt to discretize the 2d spectrum when the central abundance is known'''
    simple = np.zeros([10, 10])
    for i in range(10):
        for j in range(10):
            xmark_h = (i + 0.5) * covest[0]
            xmark_l = (i - 0.5) * covest[0]
            ymark_h = (j + 0.5) * covest[1]
            ymark_l = (j - 0.5) * covest[1]
            x_h = np.max([np.where(xax < xmark_h)])
            try:
                x_l = np.min([np.where(xax > xmark_l)])
            except ValueError:
                x_l = 0
            y_h = np.max([np.where(yax < ymark_h)])
            try:
                y_l = np.min([np.where(yax > ymark_l)])
            except ValueError:
                y_l = 0
            print(i, j, "(", x_l, ":", x_h, ") (", y_l, ":", y_h, ")",
                  xax[x_l], ":", xax[x_h], "\t", yax[y_l], ":", yax[y_h])
            simple[i][j] = np.sum(a[x_l:x_h, y_l:y_h])
    print(simple)
    return simple


if __name__ == '__main__':
    usage = "usage: %prog <input 2d matrix > -o <output file>"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("infile", type=argparse.FileType, help="input 2d matrix") 
    parser.add_argument("-o", "--output", dest="outfile", default=None,
                      help="Output file.")
    parser.add_argument("-i", "--interactive", dest="interact", default=False,
                      action="store_true", help="Interactive")
    parser.add_argument("-g", "--graphtype", dest="graphtype", default="1000",
                      help="graph type : 1000, 100, raw")
    parser.add_argument("-m", "--m", dest="m", default=1,
                      help="axis 1 decimator")
    parser.add_argument("-n", "--n", dest="n", default=1,
                      help="axis 2 decimator")
    parser.add_argument("-x", "--xlabel", dest="xlabel",
                      default=None, help="xaxis label")
    parser.add_argument("-y", "--ylabel", dest="ylabel",
                      default=None, help="yaxis label")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                      default=True, help="Verbose [default off]")

    args = parser.parse_args()
    infile = args.infile
    if not (infile and os.path.isfile(infile)):
        parser.error("Missing input file")
    if opts.xlabel != None:
        xlabel = opts.xlabel
    else:
        assert 0   # X label -x is mandatory
    if opts.ylabel != None:
        ylabel = opts.ylabel
    else:
        assert 0   # Y label -y is mandatory
    outfile = opts.outfile
    if opts.outfile == None:
        outfile = infile

    if opts.verbose:
        sys.stdout.write("Trying to read %s with numpy...\n" % infile)
    data = np.loadtxt(infile)
    print("data shape ", data.shape)
    xax = data[0, :]
    yax = data[:, 0]
    d = data[1:, 1:]
    print("d shape", d.shape)
    m = int(opts.m)
    n = int(opts.n)
    if m == 0:
        m = max(1, int(d.shape[0] / 100))
    if n == 0:
        n = max(1, int(d.shape[1] / 100))
    reduceddata = decimate(d, [m, n])
    reducedata = 1
    np.savetxt(infile + ".all.csv", reduceddata)
    print("Shape: ", reduceddata.shape)
    plt.imshow(np.log(reduceddata.T), interpolation="nearest",
               aspect="auto", origin="lower")
    plt.title(infile, fontsize=18)
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    if opts.m != "0":
        plt.xlim([0, 1000])
    if opts.n != "0":
        plt.ylim([0, 1000])
    if opts.interact:
        plt.show()
    plt.savefig(outfile + ".al1.png")

    if 0:
        sim = decimate_withcovest(d, (84.5, 80.3))
        plt.clf()
        plt.imshow(np.log(sim), interpolation="nearest",
                   aspect="auto", origin="lower")
        plt.show()
        np.savetxt("output.csv", sim)

#!/usr/bin/env python
'''This is a curve-fitting tool for interpreting kmer spectra'''
import numpy as np
import sys, os
from scipy import stats
from scipy.optimize import leastsq
from optparse import OptionParser
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from ksatools import pad

def weightedleastsquares(parameters, yvalues, xvalues, order=None):
    '''Returns vector of weighted residuals.  Used for initial fits.'''
    err = (yvalues-fitfn(xvalues, parameters, order)) / np.sqrt(yvalues+1)
    if OPTS.constrained:
        err = err + np.ones(err.shape) * sum(parameters[2:] *parameters[2:] * (parameters[2:] < 0))
    return err

def loglikelihood(parameters, yvalues, xvalues, order=None):
    '''Calculate poisson likelihood given the data and a call to fitfn.
    Returns a vector of quasi-residuals'''
    err2 = stats.poisson.logpmf(yvalues, fitfn(xvalues, parameters, order))
    for i in range(0, len(err2)):
        if err2[i] < -1000:
            err2[i] = -1000
        if np.isnan(err2[i]) or np.isinf(err2[i]):
            err2[i] = -1001
    sqerr2 = np.sqrt(-err2)
    if OPTS.verbose:
        sys.stderr.write("ll called %.1f \n"% sum(-err2))
    if OPTS.constrained:
        sqerr2 = sqerr2 + np.ones(sqerr2.shape) * sum(parameters[2:] *parameters[2:] * (parameters[2:] < 0))
    return sqerr2

def sumscorefn(parameters):
    '''Return single-value sum of a SCOREFN evaluation.  Used for derivative.'''
    xvalues = x   # global access
    yvalues = y   # global access
    return sum(SCOREFN(parameters, yvalues, xvalues, order=fittermsorder)**2)

def calculatederivative(vector, fun):
    '''Evalutes numerical second derivative of single-valued
    function (sumscorefn).'''
    import numdifftools as nd
    Hfun = nd.Hessian(fun)
    Hess = Hfun(vector)
#    plt.imshow(np.log(abs(Hess)), interpolation="nearest"); plt.colorbar(); plt.show()
    print "Hessian", Hess
    print "11", 1/np.sqrt(Hess[0][0])
    print "22", 1/np.sqrt(Hess[1][1])
    print "33", 1/np.sqrt(Hess[2][2])
    print "frac11", 1/np.sqrt(Hess[0][0])/vector[0]
    print "frac22", 1/np.sqrt(Hess[1][1])/vector[1]
    print "frac33", 1/np.sqrt(Hess[2][2])/vector[2]
    return Hess

def dumpparameters(parameters):
    '''Prints parameters.'''
    for k2 in range(len(parameters)):
        print "%d\t%.4f " % (k2, parameters[k2])
    return

def nbinompdf(xvalues, poissonlambda, alpha):
    '''Negative binomial with lambda, alpha parameterization.'''
    return np.exp(
                  stats.nbinom.logpmf(xvalues, 1/alpha, 1 / (1 + poissonlambda * alpha)))

def pevaln(xvalues, parameters, multipleorder=None):
    '''Main fit function, returns vector.'''
    if multipleorder == None:
        multipleorder = range(1, len(parameters) - 2 + 1)
    cov = np.max([0, parameters[0]])
    shap = np.max([0, parameters[1]])
    total = 0
    for i in range(0, len(multipleorder)):
        n = multipleorder[i]
        total = total + np.max([parameters[i+2], 0]) * \
            nbinompdf(xvalues, n * cov, shap/n)
    return total

def windowmask(xxx, yyy, covest, multipleorder):
    '''returns two vectors masked by specified limits.'''
    lim1 = multipleorder[0] * .5 * covest
    lim2 = (multipleorder[-1] + multipleorder[0] * .5) * covest
    lim1 = float(lim1)
    lim2 = float(lim2)
    assert lim1 >= 1 and lim2 >= 3, "Fit has run away, refusing to mask almost all the data"
    xvalues = np.array(xxx)
    temp = np.zeros(xvalues.shape)
    for i in range(len(multipleorder)):
        if multipleorder[i] == 1:
            window = .5
        else:
            window = 0.5
        lim1 = (-window + float(multipleorder[i])) * covest
        lim2 = (+window + float(multipleorder[i])) * covest
        temp = (xvalues >= max(lim1, LOWCUTOFF)) * (xvalues <= lim2) + temp
    index = np.where(temp > 0)
    returnx = xxx[index]
    returny = yyy[index]
    if len(returnx) == 0:  # There is no data in range.   This is not going to work.
        returnx = np.arange(int(lim1), int(lim2)+1)
        returny = np.zeros(returnx.shape)
    print "windowmask: returnx size", len(returnx), "max", max(returnx), "min", min(returnx)
    return returnx, returny

def writefile():
    outf = open("%s.fit.csv" % os.path.basename(OUTFILE), "w")
    outf.write("file\t%s\ncmd\t%s\ncov\t%.1f\nshape\t%.2f\ngsize\t%.1f\nngthalf\t%d\n"%
               (INFILE, " ".join(sys.argv[:]), COVERAGE, SHAPE, GENOMESIZE, NUMGTHALF))
    for k3 in range(2, len(plsq)):
        outf.write("%d\t%.1f\n" % (fittermsorder[k3-2], np.max([plsq[k3], 0])))
    outf.write("sumerr\t%f\n" % sumscorefn(plsq))
    outf.close()

def writedetails():
    OUTPUTMATRIX = np.hstack((dispx.reshape((len(dispx), 1)),
                              dispy.reshape((len(dispy), 1)),
                              model.reshape((len(model), 1))))
    np.savetxt("%s.fit.detail.csv" % OUTFILE, OUTPUTMATRIX,
               fmt=['%d', '%.1f', '%.1f'], delimiter="\t")
    print "sumerr\t%f" % sumscorefn(plsq)

def plotfit():
    '''Produce plot of fit'''
    partialx = x
    partialy = y * x
    partialm = x * fitfn(x, plsq, fittermsorder)
    plt.loglog(dispx, dispy, ".b", label="data")
    plt.loglog(dispx, model, '.-g', label="neg binom fit")
    plt.loglog(partialx, partialm, '.r', label='neg binom Fit')
    BOT = 10** int(np.log(min(x)))
    TOP = 10** int(np.log(max(x)+1))
    plt.ylim((10, 10**7))
    plt.legend(loc="lower left")
    plt.xlabel('kmer coverage')
    plt.ylabel('Number of reads ')
    plt.title('fit to %s' % (INFILE,))
    plt.text(BOT*10, 1E6, "%s = %.2f Mbases" % ("genomesize", GENOMESIZE / 1E6))
    plt.text(BOT*10, 6E5, "%s = %.1f x" % ("coverage      ", COVERAGE))
    plt.savefig("%s.fit.png" % os.path.basename(OUTFILE))
    if(OPTS.interactive):
        plt.show()

if __name__ == '__main__':
    USAGE = "usage: kmerspectrumanalyzer.py [options] <input table filename> "
    PARSER = OptionParser(USAGE)
    PARSER.add_option("-q", "--likelihood", dest="likelihood",
                      action="store_true", default=False,
                      help="Use likeliehood (slower)")
    PARSER.add_option("-g", "--guess", dest="guess",
                      default=None, help="Initial coverage guess (overrides auto-guessing)")
    PARSER.add_option("-l", "--lowcutoff", dest="lowcutoff",
                      default=10, help="Low-coverage soft cutoff (default 10)")
    PARSER.add_option("-e", "--errorbars", dest="errorbars",
                      action="store_true", default=False, help="Estimate uncertainty ")
    PARSER.add_option("-n", "--num", dest="num",
                      default=10, help="number of multiplicity terms to fit")
    PARSER.add_option("-b", "--bypass", dest="bypass",
                      default=2, help="second peak to fit for manual reordering")
    PARSER.add_option("-c", "--cutoff", dest="cutoff",
                      default="0,0", help="max,min manual hard coverage cutoffs")
    PARSER.add_option("-i", "--interactive", dest="interactive",
                      action="store_true", default=False, help="interactive plot")
    PARSER.add_option("-p", "--positiveconstrained", dest="constrained",
                      action="store_false", default=True, help="constrained fit")
    PARSER.add_option("-o", "--outstem", dest="outstem",
                      default=None, help="output file stem")
    PARSER.add_option("-v", "--verbose", dest="verbose", action="store_true",
                      default=False, help="verbose")

    (OPTS, ARGS) = PARSER.parse_args()
    try:
        INFILE = ARGS[0]
    except IndexError:
        PARSER.error("Missing table filename input argument \n%s\n"%USAGE)
    if not (INFILE and os.path.isfile(INFILE)):
        PARSER.error("Missing input filename\n%s\n"%USAGE)
    LOWCUTOFF = int(OPTS.lowcutoff)
    MAXFEV = 2000
    DATA = np.loadtxt(INFILE) # loaddata
    ORIGY = DATA[:, 1]
    ORIGX = DATA[:, 0]
    NUMBEROFTERMS = int(OPTS.num)
    SECONDPEAK = int(OPTS.bypass)
    OUTFILE = OPTS.outstem
    if OUTFILE == None:
        OUTFILE = INFILE
    print "INFILE", INFILE
    print "OUTFILE", OUTFILE
    MAX_COV, MIN_COV = [int(cov) for cov in OPTS.cutoff.split(",")]
#  pad data with zeroes representing bins with no counts.
    print "Padding... originally %d, max %d " % (len(ORIGX), max(ORIGX))
    (paddedx1, paddedy1) = pad(ORIGX.tolist(), ORIGY.tolist())
    print "Padding... finally %d "% len(paddedx1)
    sys.stdout.flush()
    padx = np.array(paddedx1)
    pady = np.array(paddedy1)
    padp = padx * pady
    lenpad = len(padp)
#  guess method 1:  mean + total of the truncated distribution
    z1 = padp[(LOWCUTOFF-1):]
    guessx = sum(padp[LOWCUTOFF-1:]) / sum(paddedy1[LOWCUTOFF-1:])
    guessy = sum(paddedy1[LOWCUTOFF-1:])
    print "Guessing... cov=", guessx, " size=", guessy

#  guess method 2: find maximum of truncated data
    maxindex = int(np.mean(np.where(padp[LOWCUTOFF-1:] == padp[LOWCUTOFF-1:].max())[0]))
    if maxindex < LOWCUTOFF - 1:
        print "WARNING! maximum %d is lower than lowcutoff %d" % (maxindex, LOWCUTOFF)
    guessmax = padx[maxindex]
    sys.stdout.write("guessmax, %d\n" % guessmax)
# guess method 3:  override guess if initial guess provided
    if OPTS.guess:
        overrideguess = int(OPTS.guess)

# Use initial guess to mask data in window before preliminary fit.
    x, y = windowmask(padx, pady, guessx, (1, 1.5))
    parms0 = np.array([guessx, .01, guessy])

#    preliminary fit
    fitfn = pevaln         # functional form of the model
    SCOREFN = weightedleastsquares

    plsq0 = leastsq(SCOREFN, parms0, args=(y, x, [1]), maxfev=MAXFEV)[0]
    print "First fit "
    dumpparameters(plsq0)

    if OPTS.likelihood:
        SCOREFN = loglikelihood
        plsq0 = leastsq(SCOREFN, plsq0, args=(y, x, [1]), maxfev=MAXFEV)[0]
        print "Second fit "
        dumpparameters(plsq0)
    else:
        SCOREFN = weightedleastsquares

    covguess = plsq0[0]
    plsq = plsq0  # initialize results in case the loop does not run
    parms0 = plsq0
# Construct list of terms order
    fittermsorder = range(1, NUMBEROFTERMS+1)
    try:
        ind = fittermsorder.index(SECONDPEAK)
        del(fittermsorder[ind])
    except ValueError:
        if len(fittermsorder) != 1:  # don't ever delete the entire list!
            del(fittermsorder[-1])
    fitterms = [fittermsorder[0], SECONDPEAK]
    fitterms.extend(fittermsorder[1:])
    fittermsorder = fitterms[0:NUMBEROFTERMS]
    print fittermsorder
    print "lenfittermsorder", len(fittermsorder)

#   successively add parameters to the fit
#    plsq = plsq[0:2]  # throw away the initial least squares fit
    for k in range(0, len(fittermsorder)):
        print "I: ", k,
        fittermsorderpartial = fittermsorder[0:k+1]
        print fittermsorderpartial
        print "applying windowmask for fit with guess", covguess, fittermsorderpartial
        x, y = windowmask(padx, pady, covguess, fittermsorderpartial)
        if OPTS.verbose:
            print "size parms0= %d  size x, y = %d" % (len(parms0), len(x))
        if k != 0:   # Do not add the next term the first time
            print "masking for guess", covguess, [fittermsorder[k]]
            tempx, tempy = windowmask(padx, pady, covguess, [fittermsorder[k]])
            tempguess = np.sum(tempy)
            print "  resulting guess ", tempguess
            parms0 = np.concatenate([plsq, [tempguess]])
        print "committing search on", fittermsorderpartial
        plsq = leastsq(SCOREFN, parms0, args=(y, x, fittermsorderpartial), diag=(1/ np.sqrt(np.minimum(np.maximum(np.abs(parms0), 10), 1E7))), maxfev=MAXFEV)[0]
        if OPTS.verbose:
            sys.stderr.write("fitting with %d params\n"%(len(parms0)))
        print "fit %.1f, %.1f, %.3f " % (plsq[0], plsq[1], plsq[2])
        covguess = plsq[0]
        print "Final fit parameters:"
        dumpparameters(plsq)
        print "*************"
    if OPTS.errorbars:
        hessian = calculatederivative(plsq, sumscorefn)
#  Display summary of the fit
    COVERAGE = plsq[0]
    SHAPE = plsq[1]
    GENOMESIZE = 0
    NUMGTHALF = np.sum(pady[np.where(padx > plsq[0] * 0.5)])
    print "parameters"
    for k in range(2, len(plsq)):
        print "%d\t%.1f" %       (fittermsorder[k-2], np.max([plsq[k], 0]))
        GENOMESIZE = GENOMESIZE + fittermsorder[k-2] * np.max([plsq[k], 0])
    print "Filename: %s"% INFILE
    print "Genomesize %.1f " % GENOMESIZE
    print "coverage   %.1f " % COVERAGE
    print "shape      %.2f " % SHAPE
    print "numdistinct %.1f " % (np.sum(np.maximum(plsq[2:], 0)))
    print "numgthalf   %d " % (NUMGTHALF)

# plot data + model
    dispx = ORIGX                          # complete x
    dispy = ORIGX * ORIGY                  # complete y data
    model = dispx * fitfn(dispx, plsq, fittermsorder)    # complete model
    writefile()
    writedetails()
    plotfit()

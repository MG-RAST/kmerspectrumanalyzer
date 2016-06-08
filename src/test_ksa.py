
from nose.tools import assert_equal, assert_almost_equal, assert_true, \
    assert_false, assert_raises, assert_is_instance

from subprocess import call

from ksatools import calccumsum, loadfile, renyispectrum, pad, smoothspectrum, calcmedian, cleanlabel, getmgrkmerspectrum, printstats, printstratify, makegraphs

import numpy as np
from numpy import array

FIXTURE1 = np.array([[2000, 1000000],
                     [20000, 10000],
                     [200000, 100],
                     [1000000, 1],
                     [2000000, 1]], dtype="float")
def test_calccumsum_pass():
    data = FIXTURE1
    cn, c1, yd, yo, zd, zo = calccumsum(data)
    pass

def test_calccumsum_00_file():
    data = loadfile("../test/test00.21")
    cn, c1, yd, yo, zd, zo = calccumsum(data)
    assert_equal(cn, 1)
    assert_equal(c1, 1)
    assert_equal(yd, 1)
    assert_equal(yo, 1)
    assert_equal(zd, 1)
    assert_equal(zo, 1)

def test_calccumsum_01_array_int():
    data = np.array([[1, 910], [2, 45]], dtype="int")
    cn, c1, yd, yo, zd, zo = calccumsum(data)
    assert_true(np.all(cn == array([1, 2])))
    assert_true(np.all(c1 == array([910, 45])))
    assert_true(np.all(yd == array([955, 45])))
    assert_true(np.all(yo == array([1000, 90])))
    assert_true(np.all(zd == array([910, 955])))
    assert_true(np.all(zo == array([910, 1000])))

def test_calccumsum_01_array_float():
    data = np.array([[1, 910], [2, 45]], dtype="float")
    cn, c1, yd, yo, zd, zo = calccumsum(data)
    assert_true(np.all(cn == array([1, 2])))
    assert_true(np.all(c1 == array([910, 45])))
    assert_true(np.all(yd == array([955, 45])))
    assert_true(np.all(yo == array([1000, 90])))
    assert_true(np.all(zd == array([910, 955])))
    assert_true(np.all(zo == array([910, 1000])))

def test_calccumsum_int():
    data = FIXTURE1
    cn, c1, yd, yo, zd, zo = calccumsum(data)
    pass

def test_calccumsum_float():
    data = FIXTURE1
    cn, c1, yd, yo, zd, zo = calccumsum(data)
    pass

def test_reyni_returns():
    exponents = np.arange(.01, 5, .01)
    renyi_sepctrum = renyispectrum(exponents, FIXTURE1)

def test_pad():
    xpadded, ypadded = pad(FIXTURE1[:,0], FIXTURE1[:,1])
    assert_true(len(xpadded) > max(FIXTURE1[:,0]) * 1.2)
    assert_true(len(ypadded) > max(FIXTURE1[:,0]) * 1.2)
    assert_true(sum(ypadded) == sum(FIXTURE1[:,1]))

def test_smooth():
    s = smoothspectrum(FIXTURE1)

def test_calcmedian():
    yd = np.array([1,2,3,4])
    y = np.array([10,20,30,40])
    num = 35
    s = calcmedian( yd, y, num )
    print s
    assert_true(s == 3.5)

def test_cleanlabel():
    assert_true( cleanlabel("file.histhist") == "file")
    assert_true( cleanlabel("file.fastq") == "file")
    assert_true( cleanlabel("file_info_contigstats.txt") == "file")
    assert_true( cleanlabel("file.stats.txt") == "file")
    assert_true( cleanlabel("file.txt") == "file")
    assert_true( cleanlabel("file.csv") == "file")
    assert_true( cleanlabel("file.037.kmerhistogram") == "file")

def test_getmgrkmerspectrum():
    g = getmgrkmerspectrum("mgm4441680.3")
    assert_true( g.shape == (55,2))

def test_getmgrkmerspectrum_fail():
    g = getmgrkmerspectrum("mgm4441680.3", mgrkey="INVALID")

def test_printstats():
    p = printstats(FIXTURE1, "out.txt")

def test_printstratify():
    p = printstratify(FIXTURE1)

def test_makesgraphs():
    for i in range(-3,27):
        print "testing visualization {:d}".format(i)
        p = makegraphs(FIXTURE1, "filename", option = i)

def test_makegraphs_d():
    p = makegraphs(FIXTURE1, "filename", option=6, dump=True)

def test_makegraphs_l():
    p = makegraphs(FIXTURE1, "filename", option=6, label="ridiculouslabel")

def test_makegraphs_30():
    p = makegraphs(FIXTURE1, "filename", option=30)


def test_cmdline_1():
    call(["plotkmerspectrum.py", "-g", "-1", "../test00.21"])

def test_cmdline_2():
    call(["plotkmerspectrum.py", "-g", "0", "../test00.21"])

#def test_cmdline_getmgr():
#    call(['plotkmerspectrum.py','-l','../test/mgrlist','-i','-g','6','-t','mgm'])

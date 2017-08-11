
from nose.tools import assert_equal, assert_almost_equal, assert_true, \
    assert_false, assert_raises, assert_is_instance

from subprocess import call

from ksatools.ksatools import loadfile, calcmedian, calccumsum, TESTDIR

import numpy as np
from numpy import array


def test_01():
    data = loadfile(TESTDIR + "test01.21")
    cn, c1, yd, yo, zd, zo = calccumsum(data)

    pass

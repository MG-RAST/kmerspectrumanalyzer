
from nose.tools import assert_equal, assert_almost_equal, assert_true, \
    assert_false, assert_raises, assert_is_instance

from ksatools import calccumsum, loadfile

import numpy as np
from numpy import array

def test_calccumsum_pass():
    data = np.array( [[2000, 1000000] , [20000, 10000], [200000, 100], [1000000, 1], [2000000, 1] ], dtype="float")
    cn, c1, yd, yo, zd, zo = calccumsum(data)
    pass

def test_calccumsum_00_file():
    data = loadfile("../test/test00.21") 
    cn, c1, yd, yo, zd, zo = calccumsum(data)
    exp = 0
    assert_equal(cn, 1)
    assert_equal(c1, 1)
    assert_equal(yd, 1)
    assert_equal(yo, 1)
    assert_equal(zd, 1)
    assert_equal(zo, 1)

def test_calccumsum_01_array_int():
    data = np.array( [ [1,910], [2,45] ], dtype="int")
    cn, c1, yd, yo, zd, zo = calccumsum(data)
    assert_true(np.all(cn ==  array([1, 2])) )
    assert_true(np.all(c1 == array([910,  45])))
    assert_true(np.all(yd == array([955,  45])))
    assert_true(np.all(yo == array([1000,  90])))
    assert_true(np.all(zd == array([910, 955]) ))
    assert_true(np.all(zo == array([910, 1000])))

def test_calccumsum_01_array_float():
    data = np.array( [ [1,910], [2,45] ], dtype="float")
    cn, c1, yd, yo, zd, zo = calccumsum(data)
    assert_true(np.all(cn ==  array([1, 2])) )
    assert_true(np.all(c1 ==  array([910,  45])))
    assert_true(np.all(yd == array([955,  45])))
    assert_true(np.all(yo == array([1000,  90])))
    assert_true(np.all(zd == array([910, 955]) ))
    assert_true(np.all(zo == array([910, 1000])))

def test_calccumsum_int():
    data = np.array( [[2000, 1000000] , [20000, 10000], [200000, 100], [1000000, 1], [2000000, 1] ], dtype="int")
    cn, c1, yd, yo, zd, zo = calccumsum(data)
    pass

def test_calccumsum_float():
    data = np.array( [[2000, 1000000] , [20000, 10000], [200000, 100], [1000000, 1], [2000000, 1] ], dtype="float")
    cn, c1, yd, yo, zd, zo = calccumsum(data)
    pass

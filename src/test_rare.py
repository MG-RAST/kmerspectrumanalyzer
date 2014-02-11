
from nose.tools import assert_equal, assert_almost_equal, assert_true, \
    assert_false, assert_raises, assert_is_instance
import numpy as np
from numpy import array

from rare import fract

# fract(aa, epsilon, threshold)
def test_calccumsum_one_full():
    data = np.array( [[100, 1000] ], dtype="float" )
    result = fract(data, 1, 1) 
    assert_equal(result , 1)

def test_calccumsum_one_zero():
    data = np.array( [[100, 1000] ], dtype="float" )
    result = fract(data, 1, 100000) 
    assert_equal(result , 0)

def test_calccumsum_one_half():
    data = np.array( [[100, 1000] ], dtype="float" )
    result = fract(data, 1, 100) 
    assert_true(result > 0.3)
    assert_true(result < 0.7)

def test_calccumsum_two_full():
    data = np.array( [[100, 1000] , [20000, 10000] ], dtype="float" )
    result = fract(data, 1, 1) 
    assert_equal(result, 1)
 
def test_calccumsum_two_zero():
    data = np.array( [[100, 1000] , [20000, 10000] ], dtype="float" )
    result = fract(data, 1, 100000) 
    assert_equal(result, 0)

def test_calccumsum_two_tenth():
    data = np.array( [[100, 1000] , [1000, 1000] ], dtype="float" )
    result = fract(data, 1, 200) 
    assert_true(result > 0.88)
    assert_true(result < 0.92) 


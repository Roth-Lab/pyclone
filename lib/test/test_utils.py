'''
Created on 2012-09-17

@author: andrew
'''
from __future__ import division

import unittest

import pyclone.utils as utils

from math import log, exp

class Test(unittest.TestCase):
    def test_log_sum_exp(self):
        X = [100, 200, 300]
        
        log_X = [log(x) for x in X]
        
        s_X = sum(X)
        
        log_s_X = utils.log_sum_exp(log_X)
        
        self.assertEqual(exp(log_s_X), s_X)
    
    def test_log_space_normalise(self):
        X = [100, 200, 300]
        
        log_X = [log(x) for x in X]
        
        norm_X = [x / sum(X) for x in X]
        
        log_norm_X = utils.log_space_normalise(log_X)
        
        self.assertAlmostEqual([exp(x) for x in log_norm_X][0], norm_X[0], 6)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()

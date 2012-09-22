'''
Created on 2012-09-19

@author: andrew
'''
from __future__ import division

import unittest

from pyclone.utils import binomial_rvs, discrete_rvs

class Test(unittest.TestCase):
    def test_binomial(self):
        n = 10
        p = 0.001
        
        x = [binomial_rvs(n, p) for _ in range(1000000)]
        
        mean = sum(x) / len(x)
        
        self.assertAlmostEqual(mean, p * n, 2)
    
    def test_discrete(self):
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_binomial']
    unittest.main()

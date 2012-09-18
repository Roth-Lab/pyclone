'''
Created on 2012-08-22

@author: andrew
'''
from __future__ import division

import unittest

from math import log, exp

from pyclone.utils import SimpsonsRuleIntegrator

class Test(unittest.TestCase):
    def test_simple_integral(self):
        a = 1
        b = 2
        
        integrator = SimpsonsRuleIntegrator(a, b, 1000)
        
        # log(x**5)
        def log_f(x):
            if x == 0:
                return float('-inf')
            else:
                return 15 * log(x)

        observed_value = exp(integrator.log_integrate(log_f))
        
        F = lambda x: 1 / 16 * x ** 16
        
        expected_value = F(b) - F(a)
        
        self.assertAlmostEqual(observed_value, expected_value, delta=1e-6)

    def test_odd_mesh_size_raises_error(self):
        self.assertRaises(Exception, SimpsonsRuleIntegrator, 0, 1, 99)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()

'''
Created on 2012-08-22

@author: andrew
'''
from __future__ import division

import unittest

from math import log, exp

from pyclone._likelihoods import BetaBinomialLikelihood, DataPoint
from pyclone.utils import SimpsonsRuleIntegrator, log_integrate


class Test(unittest.TestCase):
    def test_simple_integral(self):
        a = 0
        b = 2
        
        integrator = SimpsonsRuleIntegrator(a, b, 100)
        
        # log(x**2)
        def log_f(x):
            if x == 0:
                return float('-inf')
            else:
                return 5 * log(x)

        observed_value = exp(integrator.log_integrate(log_f))
        
        F = lambda x: 1 / 6 * x ** 6
        
        expected_value = F(b) - F(a)
        
        self.assertAlmostEqual(observed_value, expected_value, delta=1e-6)

        print exp(log_integrate(log_f, a, b, 10)), expected_value
        

    def test_odd_mesh_size_raises_error(self):
        self.assertRaises(Exception, SimpsonsRuleIntegrator, 0, 1, 99)
    
    def test_likelihood(self):
        phi = 0.01
        
        data_point = DataPoint(1000, 2000, [0.99, ], [0.01, 0.5], [1, ], [10, 10])        
        
        likelihood = BetaBinomialLikelihood(data_point, mesh_size=10)
        
        print likelihood.compute_log_likelihood(phi)
        
        likelihood = BetaBinomialLikelihood(data_point, mesh_size=100)
        
        print likelihood.compute_log_likelihood(phi)

        likelihood = BetaBinomialLikelihood(data_point, mesh_size=500)
        
        print likelihood.compute_log_likelihood(phi)

        likelihood = BetaBinomialLikelihood(data_point, mesh_size=1000)
        
        print likelihood.compute_log_likelihood(phi)                

        likelihood = BetaBinomialLikelihood(data_point, mesh_size=10000)
        
        print likelihood.compute_log_likelihood(phi)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()

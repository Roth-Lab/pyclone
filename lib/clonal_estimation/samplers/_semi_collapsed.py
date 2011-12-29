'''
Created on 2011-12-29

@author: Andrew Roth
'''
from __future__ import division

from math import log
from numpy.random import binomial, beta
from random import random

from clonal_estimation.model import SemiCollapsedLikelihood, DataPoint
from clonal_estimation.utils import log_binomial_pdf

class SemiCollapsedSampler(object):
    def run(self, data_point, max_iters=100000, thin=1, burnin=0):
        results = {'phi' : [], 'd_v' : []}
        
        likelihood = SemiCollapsedLikelihood(data_point)
        
        phi = random()
        
        d_v = binomial(data_point.d, phi)
        
        for i in range(max_iters):
            d_v = self._update_d_v(likelihood, d_v, phi)
            phi = self._update_phi(likelihood, d_v, phi)
            
            if i % thin == 0 and i >= burnin:
                results['phi'].append(phi)
                results['d_v'].append(d_v)
        
                print i, d_v, phi
        
        return results
    
    def _update_d_v(self, likelihood, old_d_v, phi):
        d = likelihood.d
        
        new_d_v = binomial(d, phi)
        
        numerator = likelihood.compute_log_likelihood(new_d_v, phi) + log_binomial_pdf(old_d_v, d, phi)
        
        denominator = likelihood.compute_log_likelihood(old_d_v, phi) + log_binomial_pdf(new_d_v, d, phi)
        
        log_ratio = numerator - denominator
        
        u = random()
        
        if log_ratio >= log(u):
            return new_d_v
        else:
            return old_d_v
    
    def _update_phi(self, likelihood, d_v, old_phi):
        d = likelihood.d
        
        d_r = d - d_v
        
        a = d_v + 1
        b = d_r + 1
        
        new_phi = beta(a, b)
        
        return new_phi 

'''
Created on 2011-12-29

@author: Andrew Roth
'''
from __future__ import division

from collections import OrderedDict
from math import log

from pyclone.utils import log_sum_exp, log_binomial_coefficient, log_binomial_likelihood

class DataPoint(object):
    def __init__(self, a, d, mu_r, mu_v, log_pi_r, log_pi_v):
        self.a = a
        self.d = d
        
        self.log_pi_r = log_pi_r
        self.log_pi_v = log_pi_v
        
        self.mu_r = mu_r
        self.mu_v = mu_v

#=======================================================================================================================
# Likelihood
#=======================================================================================================================
class BinomialLikelihood(object):
    def __init__(self, data_point):
        self.a = data_point.a
        self.d = data_point.d
        
        self.mu_r = data_point.mu_r
        self.mu_v = data_point.mu_v
        
        self.log_pi_r = data_point.log_pi_r
        self.log_pi_v = data_point.log_pi_v
        
        self._ll_cache = OrderedDict()
        
        self._log_binomial_norm_const = log_binomial_coefficient(self.d, self.a)
        
    def compute_log_likelihood(self, phi):        
        if phi not in self._ll_cache:
            self._ll_cache[phi] = self._log_likelihood(phi)
            
            if len(self._ll_cache) > 1000:
                self._ll_cache.popitem(last=False)
        
        return self._ll_cache[phi]
    
    def _log_likelihood(self, phi):    
        ll = []
        
        for mu_r, log_pi_r in zip(self.mu_r, self.log_pi_r):
            for mu_v, log_pi_v in zip(self.mu_v, self.log_pi_v):
                temp = log_pi_r + log_pi_v + self._log_complete_likelihood(phi, mu_r, mu_v)
                
                ll.append(temp)
        
        return log_sum_exp(ll)
    
    def _log_complete_likelihood(self, phi, mu_r, mu_v):
        mu = (1 - phi) * mu_r + phi * mu_v
        
        return self._log_binomial_norm_const + log_binomial_likelihood(self.a, self.d, mu)

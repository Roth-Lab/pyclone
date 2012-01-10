'''
Created on 2011-12-29

@author: Andrew Roth
'''
from __future__ import division

from collections import OrderedDict
from math import log

from pyclone.utils import log_sum_exp, log_binomial_coefficient, log_binomial_likelihood

class DataPoint(object):
    def __init__(self, a, d, mu_n, mu_r, mu_v, pi_r, pi_v):
        self.a = a
        self.d = d
        
        self.pi_r = pi_r
        self.pi_v = pi_v
        
        self.mu_n = mu_n
        self.mu_r = mu_r
        self.mu_v = mu_v

#=======================================================================================================================
# Likelihood
#=======================================================================================================================
class BinomialLikelihood(object):
    def __init__(self, data_point):
        self.a = data_point.a
        self.d = data_point.d
        
        self.mu_n = data_point.mu_n
        self.mu_r = data_point.mu_r
        self.mu_v = data_point.mu_v
                
        self.pi_r = data_point.pi_r
        self.pi_v = data_point.pi_v
        
        self._ll_cache = OrderedDict()
        
        self._log_binomial_norm_const = log_binomial_coefficient(self.d, self.a)
        
    def compute_log_likelihood(self, phi, s):
        key = (phi, s)
        
        if key not in self._ll_cache:
            self._ll_cache[key] = self._log_likelihood(phi, s)
            
            if len(self._ll_cache) > 1000:
                self._ll_cache.popitem(last=False)
        
        return self._ll_cache[key]
    
    def _log_likelihood(self, phi, s):    
        ll = []
        
        for mu_r, pi_r in zip(self.mu_r, self.pi_r):
            for mu_v, pi_v in zip(self.mu_v, self.pi_v):
                temp = log(pi_r) + log(pi_v) + self._log_complete_likelihood(phi, s, mu_r, mu_v)
                
                ll.append(temp)
        
        return log_sum_exp(ll)
    
    def _log_complete_likelihood(self, phi, s, mu_r, mu_v):
        mu_t = (1 - phi) * mu_r + phi * mu_v
        
        mu = (1 - s) * self.mu_n + s * mu_t
        
        return self._log_binomial_norm_const + log_binomial_likelihood(self.a, self.d, mu)

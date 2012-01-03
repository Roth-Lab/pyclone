'''
Created on 2011-12-29

@author: Andrew Roth
'''
from __future__ import division

from collections import OrderedDict
from math import log

from clonal_estimation.utils import log_sum_exp, log_factorial

class DataPoint(object):
    def __init__(self, a, d, mu_r, mu_v, pi_r, pi_v):
        self.a = a
        self.d = d
        
        self.pi_r = pi_r
        self.pi_v = pi_v
        
        self.mu_r = mu_r
        self.mu_v = mu_v

class Likelihood(object):
    def __init__(self, data_point):
        self.a = data_point.a
        self.d = data_point.d
        
        self._log_mu_r = []
        self._log_mu_v = []
        
        self._log_pi_r = []
        self._log_pi_v = []
        
        for m_r, p_r in zip(data_point.mu_r, data_point.pi_r):
            if p_r == 0:
                continue
            
            print m_r, p_r
            
            self._log_pi_r.append(log(p_r))
            self._log_mu_r.append((log(m_r), log(1 - m_r)))
        
            
            
        for m_v, p_v in zip(data_point.mu_v, data_point.pi_v):
            print m_v, p_v
            
            if p_v == 0:
                continue
            
            self._log_pi_v.append(log(p_v))
            self._log_mu_v.append((log(m_v), log(1 - m_v)))
    
        self._partial_ll_cache = {}
        
        self._ll_cache = OrderedDict()
        
        self._factorial_cache = [log_factorial(x) for x in range(self.d + 1)]

    def compute_log_likelihood(self, *args):
        pass
                
    def _update_ll_cache(self, d_v):
        a = self.a
        d = self.d
        
        d_r = self.d - d_v
        
        a_min = max(0, a - d_r)
        a_max = min(a, d_v)
        
        lf = self._factorial_cache
            
        log_inner_sum = []
        
        for a_v in range(a_min, a_max + 1):     
            a_r = a - a_v
            
            temp = lf[d] - lf[a_r] - lf[d_r - a_r] - lf[a_v] - lf[d_v - a_v]
            
            temp += self._log_binomial_likelihood(a_r, d_r, self._log_mu_r, self._log_pi_r)
            
            temp += self._log_binomial_likelihood(a_v, d_v, self._log_mu_v, self._log_pi_v)
            
            log_inner_sum.append(temp)
        
        self._partial_ll_cache[d_v] = log_sum_exp(log_inner_sum)
    
    def _log_binomial_likelihood(self, a, d, log_mu, log_pi):
        temp = []
        
        for mu, pi in zip(log_mu, log_pi):
            temp.append(pi + a * mu[0] + (d - a) * mu[1])
        
        return log_sum_exp(temp)

class SemiCollapsedLikelihood(Likelihood):
    def compute_log_likelihood(self, d_v, phi):
        if d_v not in self._partial_ll_cache:
            self._update_ll_cache(d_v)
        
        d_r = self.d - d_v

        return d_v * log(phi) + d_r * log(1 - phi) + self._partial_ll_cache[d_v]    

class CollapsedLikelihood(Likelihood):
    def __init__(self, data_point):
        Likelihood.__init__(self, data_point)
        
        for d_v in range(data_point.d + 1):
            self._update_ll_cache(d_v)
            
    def compute_log_likelihood(self, phi):
        if phi in self._ll_cache:
            ll = self._ll_cache[phi]
        else:  
            temp = []
            
            log_phi = [log(phi), log(1 - phi)]
            
            for d_v in range(self.d + 1):
                d_r = self.d - d_v
                
                partial_ll = d_v * log_phi[0] + d_r * log_phi[1] + self._partial_ll_cache[d_v]
                
                temp.append(partial_ll)
            
            ll = log_sum_exp(temp)
            
            self._ll_cache[phi] = ll
            
            if len(self._ll_cache) > 1000:
                self._ll_cache.popitem(last=False)
        
        return ll

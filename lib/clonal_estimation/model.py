'''
Created on 2011-12-29

@author: Andrew Roth
'''
from collections import OrderedDict
from math import log

from clonal_estimation.utils import log_sum_exp

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
        
        self._mu_r = []
        self._mu_v = []
        
        self._log_pi_r = []
        self._log_pi_v = []
        
        for g_r, p_r in enumerate(data_point.pi_r):
            if p_r == 0:
                continue
            
            self._log_pi_r.append(log(p_r))
            self._mu_r.append(data_point.mu_r[g_r])
        
            print data_point.mu_r[g_r]
            
        for g_v, p_v in enumerate(data_point.pi_v):
            if p_v == 0:
                continue
            
            self._log_pi_v.append(log(p_v))
            self._mu_v.append(data_point.mu_v[g_v])
            
            print data_point.mu_v[g_v]
    
        self._partial_ll_cache = {}
        
        self._ll_cache = OrderedDict()

    def compute_log_likelihood(self, *args):
        pass

    def _update_ll_cache(self, d_v):
        d_r = self.d - d_v
        
        temp = []
        
        for mu_r, log_pi_r in zip(self._mu_r, self._log_pi_r):
            for mu_v, log_pi_v in zip(self._mu_v, self._log_pi_v): 
                temp.append(log_pi_r + log_pi_v + self._compute_latent_term(d_r, d_v, mu_r, mu_v))
            
        self._partial_ll_cache[d_v] = log_sum_exp(temp)    
        
    def _compute_latent_term(self, d_r, d_v, mu_r, mu_v):
        ll = []
        
        for a_r in range(d_r + 1):
            for a_v in range(d_v + 1):
                if a_r + a_v != self.a:
                    continue
                else:
                    ll.append(log_binomial_pdf(a_r, d_r, mu_r) + log_binomial_pdf(a_v, d_v, mu_v))

        return log_sum_exp(ll)

class SemiCollapsedLikelihood(Likelihood):
    def compute_log_likelihood(self, d_v, phi):
        if d_v not in self._partial_ll_cache:
            self._update_ll_cache(d_v)

        return log_binomial_pdf(d_v, self.d, phi) + self._partial_ll_cache[d_v]    

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
            
            for d_v in range(self.d + 1):
                partial_ll = self._partial_ll_cache[d_v]
                partial_ll += log_binomial_pdf(d_v, self.d, phi)
                
                temp.append(partial_ll)
            
            ll = log_sum_exp(temp)
            
            self._ll_cache[phi] = ll
            
            if len(self._ll_cache) > 1000:
                self._ll_cache.popitem(last=False)
        
        return ll
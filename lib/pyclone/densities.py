'''
Created on 2012-11-21

@author: andrew
'''
from collections import namedtuple

from pydp.densities import Density, log_binomial_pdf
from pydp.utils import log_sum_exp

CellSampledData = namedtuple('CellSampledData', ['a', 'd', 'mu_r', 'mu_v', 'log_pi_r', 'log_pi_v'])

FragmentSampledData = namedtuple('FragmentSampledData', 
                                 ['a', 'd', 'mu_r', 'mu_v', 'log_pi_r', 'log_pi_v', 'cn_r', 'cn_v'])

class CellSampledDensity(Density):
    def log_p(self, data, params):
        ll = []
        
        for mu_r, log_pi_r in zip(data.mu_r, data.log_pi_r):
            for mu_v, log_pi_v in zip(data.mu_v, data.log_pi_v):
                temp = log_pi_r + log_pi_v + self._log_binomial_likelihood(data.a, data.d, params.phi, params.s, mu_r, mu_v)
                
                ll.append(temp)
        
        return log_sum_exp(ll)
    
    def _log_binomial_likelihood(self, a, d, phi, s, mu_r, mu_v):
        mu_N = mu_r
        
        mu_T = (1 - phi) * mu_r + phi * mu_v
        
        mu = (1 - s) * mu_N + s * mu_T
        
        return log_binomial_pdf(a, d, mu)
    
class FragmentSampledDensity(Density):
    def log_p(self, data, params):
        ll = []
        
        for mu_r, log_pi_r, cn_r in zip(data.mu_r, data.log_pi_r, data.cn_r):
            for mu_v, log_pi_v, cn_v in zip(data.mu_v, data.log_pi_v, data.cn_v):
                temp = log_pi_r + log_pi_v + self._log_binomial_likelihood(data.a,
                                                                           data.d,
                                                                           params.phi,
                                                                           params.s,
                                                                           mu_r,
                                                                           mu_v,
                                                                           cn_r,
                                                                           cn_v)
                
                ll.append(temp)
        
        return log_sum_exp(ll)
    
    def _log_binomial_likelihood(self, a, d, phi, s, mu_r, mu_v, cn_r, cn_v):
        cn_N = 2
        
        # Compute probabilities of sampling from ref and var populations
        p_r = (1 - s) * cn_N + s * (1 - phi) * cn_r
        p_v = s * phi * cn_v
        
        # Normalise probabilities
        p_r = p_r / (p_r + p_v)
        p_v = 1 - p_r

        mu = p_r * mu_r + p_v * mu_v
        
        return log_binomial_pdf(a, d, mu) 
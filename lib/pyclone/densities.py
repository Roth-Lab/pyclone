'''
Created on 2012-11-21

@author: andrew
'''
from collections import namedtuple

from pydp.densities import Density, log_binomial_pdf
from pydp.utils import log_sum_exp

CellSampledData = namedtuple('CellSampledData', 
                             ['b', 'd', 'eps', 'mu_v', 'log_pi'])

FragmentSampledData = namedtuple('FragmentSampledData',
                                 ['b', 'd', 'eps', 'cn_r', 'cn_v', 'mu_v', 'log_pi'])

class CellSampledDensity(Density):
    def log_p(self, data, params):
        ll = []
        
        for mu_v, log_pi in zip(data.mu_v, data.log_pi):
            temp = log_pi + self._log_binomial_likelihood(data.b, data.d, params.phi, params.s, data.eps, mu_v)
                
            ll.append(temp)
        
        return log_sum_exp(ll)
    
    def _log_binomial_likelihood(self, a, d, phi, s, eps, mu_v):
        mu_N = eps
        
        mu_T = (1 - phi) * eps + phi * mu_v
        
        mu = (1 - s) * mu_N + s * mu_T
        
        return log_binomial_pdf(a, d, mu)
    
class FragmentSampledDensity(Density):
    def log_p(self, data, params):
        ll = []
        
        for cn_r, cn_v, mu_v, log_pi  in zip(data.cn_r, data.cn_v, data.mu_v, data.log_pi):
            temp = log_pi + self._log_binomial_likelihood(data.b,
                                                          data.d,
                                                          params.phi,
                                                          params.s,
                                                          data.eps,
                                                          mu_v,
                                                          cn_r,
                                                          cn_v)
            
            ll.append(temp)
        
        return log_sum_exp(ll)
    
    def _log_binomial_likelihood(self, b, d, phi, s, eps, mu_v, cn_r, cn_v):
        cn_N = 2
        
        # Compute probabilities of sampling from ref and var populations
        p_r = (1 - s) * cn_N + s * (1 - phi) * cn_r
        p_v = s * phi * cn_v
        
        # Normalise probabilities
        p_r = p_r / (p_r + p_v)
        p_v = 1 - p_r

        mu = p_r * eps + p_v * mu_v
        
        return log_binomial_pdf(b, d, mu) 

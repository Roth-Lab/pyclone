'''
Created on 2011-12-29

@author: Andrew Roth
'''
from __future__ import division

from math import lgamma as log_gamma

from pyclone.utils import SimpsonsRuleIntegrator, memoized, log_beta_pdf, log_binomial_likelihood, log_sum_exp

class Customer(object):
    def __init__(self):
        self.data = None
        
    

class PyCloneLikelihood(object):
    def __init__(self, a, d, ref_priors, var_priors):
        self.a = a
        self.d = d
        
        self.mu_r = ref_priors.keys()
        self.mu_v = var_priors.keys()
        
        self.log_pi_r = self._get_log_mix_weights(ref_priors.values())
        self.log_pi_v = self._get_log_mix_weights(var_priors.values())        

        
    def _get_log_mix_weights(self, delta):
        log_denominator = log_gamma(sum(delta) + 1)
        
        log_mix_weights = []
        
        for i, d_i in enumerate(delta):
            log_numerator = log_gamma(d_i + 1)
            
            for j, d_j in enumerate(delta):
                if i != j:
                    log_numerator += log_gamma(d_j)
            
            log_mix_weights.append(log_numerator - log_denominator)
        
        return log_mix_weights
    
    def _log_binomial_likelihood(self, phi, mu_r, mu_v):
        mu = (1 - phi) * mu_r + phi * mu_v
        
        return log_binomial_likelihood(self.a, self.d, mu)

class BetaBinomialLikelihood(PyCloneLikelihood):
    def __init__(self, a, d, ref_priors, var_priors):
        PyCloneLikelihood.__init__(self, a, d, ref_priors, var_priors)

        self._integrator = SimpsonsRuleIntegrator(0, 1, mesh_size=100)

    @memoized
    def evaluate(self, phi, s):
        ll = []
        
        for mu_r, log_pi_r in zip(self.mu_r, self.log_pi_r):
            for mu_v, log_pi_v in zip(self.mu_v, self.log_pi_v):
                a_v = s * mu_v
                b_v = s - a_v
                
                log_f = lambda x: self._log_binomial_likelihood(phi, mu_r, x) + log_beta_pdf(x, a_v, b_v)

                temp = log_pi_r + log_pi_v + self._integrator.log_integrate(log_f)
                
                ll.append(temp)
        
        return log_sum_exp(ll)
    
class BinomialLikelihood(PyCloneLikelihood):    
    @memoized
    def evaluate(self, phi):
        ll = []
        
        for mu_r, log_pi_r in zip(self.mu_r, self.log_pi_r):
            for mu_v, log_pi_v in zip(self.mu_v, self.log_pi_v):
                temp = log_pi_r + log_pi_v + self._log_binomial_likelihood(phi, mu_r, mu_v)
                
                ll.append(temp)
        
        return log_sum_exp(ll) 

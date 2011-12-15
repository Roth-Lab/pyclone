from __future__ import division

from math import log
from random import random

from utils import log_sum_exp, log_binomial_pdf, log_beta

def cellular_frequency_sampler(a, d, pi_r, pi_v, mu, max_iters=100000, thin=1, burnin=0):
    results = {'phi' : []}
    
    x = DataPoint(a, d, pi_r, pi_v, mu)
    
    phi = random()
    
    for i in range(max_iters):
        phi = update_phi(x, phi)
        
        if i % thin == 0 and i >= burnin:
            results['phi'].append(phi)
            print i, phi
        
    return results

def update_phi(data_point, old_phi):   
    new_phi = random()
    
    numerator = data_point.compute_log_likelihood(new_phi)
    
    denominator = data_point.compute_log_likelihood(old_phi)
    
    log_ratio = numerator - denominator
    
    u = random()
    
    if log_ratio >= log(u):
        return new_phi
    else:
        return old_phi

class DataPoint(object):
    def __init__(self, a, d, pi_r, pi_v, mu):
        self.a = a
        self.d = d
        
        self._mu_r = []
        self._mu_v = []
        
        self._log_pi_r = []
        self._log_pi_v = []
        
        for g_r, p_r in enumerate(pi_r):
            if p_r == 0:
                continue
            
            self._log_pi_r.append(log(p_r))
            self._mu_r.append(mu[g_r])
        
            print mu[g_r]
            
        for g_v, p_v in enumerate(pi_v):
            if p_v == 0:
                continue
            
            self._log_pi_v.append(log(p_v))
            self._mu_v.append(mu[g_v])
            
            print mu[g_v]
    
        self._init_ll_cache()
        
    def _init_ll_cache(self):
        self._ll_cache = {}
        
        for d_v in range(self.d + 1):
            self._ll_cache[d_v] = self._compute_partial_log_likelihood(d_v)       
    
    def compute_log_likelihood(self, phi):
        temp = []
        
        for d_v in range(self.d + 1):
            ll = self._ll_cache[d_v]
            ll += log_binomial_pdf(d_v, self.d, phi)
            
            temp.append(ll)
            
        return log_sum_exp(temp)
    
    def _compute_partial_log_likelihood(self, d_v):
        d_r = self.d - d_v
        
        temp = []
        
        for mu_r, log_pi_r in zip(self._mu_r, self._log_pi_r):
            for mu_v, log_pi_v in zip(self._mu_v, self._log_pi_v): 
                temp.append(log_pi_r + log_pi_v + self._compute_latent_term(d_r, d_v, mu_r, mu_v))
            
        return log_sum_exp(temp)    
        
    def _compute_latent_term(self, d_r, d_v, mu_r, mu_v):
        ll = []
        
        for a_r in range(d_r + 1):
            for a_v in range(d_v + 1):
                if a_r + a_v != self.a:
                    continue
                else:
                    ll.append(log_binomial_pdf(a_r, d_r, mu_r) + log_binomial_pdf(a_v, d_v, mu_v))

        return log_sum_exp(ll)
            
if __name__ == "__main__":
    mu = [0.01, 0.5, 0.99]
    pi_r = [0, 0, 1]
    pi_v = [0.2, 0.8, 0]
    
    a = 60
    d = 100
    
    results = cellular_frequency_sampler(a, d, pi_r, pi_v, mu, burnin=50000, thin=10, max_iters=100000)
    
    from math import exp
    
    import matplotlib.pyplot as plot
    
    n = 1000
    data_point = DataPoint(a, d, pi_r, pi_v, mu)
    
    x = []
    y = []
    
    for i in range(1, n):
        x.append(i / n)
        y.append(exp(data_point.compute_log_likelihood(x[-1])))
        
    
    plot.plot(x, y)
    plot.hist(results['phi'], normed=True, bins=100)
    plot.show()

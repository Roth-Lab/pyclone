from __future__ import division

from math import exp, log
from random import random

from numpy.random import binomial, beta

from utils import log_sum_exp, log_binomial_pdf

def cellular_frequency_sampler(a, d, pi_r, pi_v, mu, max_iters=100000, thin=1, burnin=0):
    results = {'phi' : [], 'd_v' : []}
    
    x = DataPoint(a, d, pi_r, pi_v, mu)
    
    phi = random()
    d_v = binomial(x.d, phi)
    
    for i in range(max_iters):
        d_v = update_d_v(x, d_v, phi)
        phi = update_phi(x, d_v, phi)
        
        if i % thin == 0 and i >= burnin:
            results['phi'].append(phi)
            results['d_v'].append(d_v)
            print i, d_v, phi
        
    return results

def update_d_v(data_point, old_d_v, phi):
    d = data_point.d
    
    new_d_v = binomial(d, phi)
    
    numerator = data_point.compute_log_likelihood(new_d_v, phi) + log_binomial_pdf(old_d_v, d, phi)
    
    denominator = data_point.compute_log_likelihood(old_d_v, phi) + log_binomial_pdf(new_d_v, d, phi)
    
    ratio = exp(numerator - denominator)
    
    u = random()
    
#    print ratio, u
    
    if ratio >= u:
        return new_d_v
    else:
        return old_d_v

def update_phi(data_point, d_v, old_phi):
    d = data_point.d
    
    d_r = d - d_v
    
    a = d_v + 1
    b = d_r + 1
    
    new_phi = beta(a, b)
    
    return new_phi

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
        
        self._ll_cache = {}
    
    def compute_log_likelihood(self, d_v, phi):
        ll = 0
        
        d_r = self.d - d_v
        
        temp = []
        
        for mu_r, log_pi_r in zip(self._mu_r, self._log_pi_r):
            for mu_v, log_pi_v in zip(self._mu_v, self._log_pi_v):
                temp.append(log_pi_r + log_pi_v + self._compute_latent_term(d_r, d_v, mu_r, mu_v))
            
        ll += log_sum_exp(temp)
        
        ll += log_binomial_pdf(d_v, self.d, phi)
        
        return ll          
        
    def _compute_latent_term(self, d_r, d_v, mu_r, mu_v):
        key = (d_r, d_v, mu_r, mu_v)
        
        if key in self._ll_cache:
            return self._ll_cache[key]
        
        ll = []
        
        for a_r in range(d_r + 1):
            for a_v in range(d_v + 1):
                if a_r + a_v != self.a:
                    ll.append(float('-inf'))
                else:
                    ll.append(log_binomial_pdf(a_r, d_r, mu_r) + log_binomial_pdf(a_v, d_v, mu_v))

        self._ll_cache[key] = log_sum_exp(ll)
        
        return self._ll_cache[key]

if __name__ == "__main__":
    mu = [0.01, 0.5, 0.99]
    pi_r = [0, 0, 1]
    pi_v = [0.2, 0.8, 0]
    
    a = 60
    d = 100
    
    results = cellular_frequency_sampler(a, d, pi_r, pi_v, mu, burnin=50000, thin=100, max_iters=1000000)
    
    from model_1_collapsed import DataPoint as CollapsedDataPoint 
    import matplotlib.pyplot as plot
    
    n = 1000
    data_point = CollapsedDataPoint(a, d, pi_r, pi_v, mu)
    
    x = []
    y = []
    
    for i in range(1, n):
        x.append(i / n)
        y.append(exp(data_point.compute_log_likelihood(x[-1])))
        
    
    plot.plot(x, y)
    plot.hist(results['phi'], normed=True, bins=100)
    plot.show()
    
    
#    x = DataPoint(a, d, pi_r, pi_v, mu)
#    
#    for i in range(1, d):
#        print x.compute_log_likelihood(1000, i / d)
#    print x.compute_log_likelihood(500, 0.99)

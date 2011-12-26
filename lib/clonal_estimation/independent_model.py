'''
Created on 2011-12-15

@author: Andrew Roth
'''
from __future__ import division

from math import log
from random import random
from collections import OrderedDict

#from numpy.random import binomial, beta

from utils import log_sum_exp, log_binomial_pdf
    
#class SemiCollapsedSampler(object):
#    def run(self, data_point, max_iters=100000, thin=1, burnin=0):
#        results = {'phi' : [], 'd_v' : []}
#        
#        likelihood = SemiCollapsedLikelihood(data_point)
#        
#        phi = random()
#        d_v = binomial(data_point.d, phi)
#        
#        for i in range(max_iters):
#            d_v = self._update_d_v(likelihood, d_v, phi)
#            phi = self._update_phi(likelihood, d_v, phi)
#            
#            if i % thin == 0 and i >= burnin:
#                results['phi'].append(phi)
#                results['d_v'].append(d_v)
#        
#                print i, d_v, phi
#        
#        return results
#    
#    def _update_d_v(self, likelihood, old_d_v, phi):
#        d = data_point.d
#        
#        new_d_v = binomial(d, phi)
#        
#        numerator = likelihood.compute_log_likelihood(new_d_v, phi) + log_binomial_pdf(old_d_v, d, phi)
#        
#        denominator = likelihood.compute_log_likelihood(old_d_v, phi) + log_binomial_pdf(new_d_v, d, phi)
#        
#        log_ratio = numerator - denominator
#        
#        u = random()
#        
#        if log_ratio >= log(u):
#            return new_d_v
#        else:
#            return old_d_v
#    
#    def _update_phi(self, likelihood, d_v, old_phi):
#        d = likelihood.d
#        
#        d_r = d - d_v
#        
#        a = d_v + 1
#        b = d_r + 1
#        
#        new_phi = beta(a, b)
#        
#        return new_phi    

class CollapsedSampler(object):
    def run(self, data_point, max_iters=100000, thin=1, burnin=0):
        results = {'phi' : []}

        likelihood = CollapsedLikelihood(data_point)
            
        phi = random()
        
        for i in range(max_iters):
            phi = self._update_phi(likelihood, phi)
            
            if i % thin == 0 and i >= burnin:
                results['phi'].append(phi)
                print i, phi
        
        return results

    def _update_phi(self, likelihood, old_phi):   
        new_phi = random()
        
        numerator = likelihood.compute_log_likelihood(new_phi)
        
        denominator = likelihood.compute_log_likelihood(old_phi)
        
        log_ratio = numerator - denominator
        
        u = random()
        
        if log_ratio >= log(u):
            return new_phi
        else:
            return old_phi

#=======================================================================================================================
# Likelihood functions.
#=======================================================================================================================
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
            self._mu_r.append(data_point.mu[g_r])
        
            print data_point.mu[g_r]
            
        for g_v, p_v in enumerate(data_point.pi_v):
            if p_v == 0:
                continue
            
            self._log_pi_v.append(log(p_v))
            self._mu_v.append(data_point.mu[g_v])
            
            print data_point.mu[g_v]
    
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

#=======================================================================================================================
# Data Point
#=======================================================================================================================
class DataPoint(object):
    def __init__(self, a, d, pi_r, pi_v, mu):
        self.a = a
        self.d = d
        
        self.pi_r = pi_r
        self.pi_v = pi_v
        
        self.mu = mu
    
            
if __name__ == "__main__":
    mu = [0.01, 0.5, 0.99]
    pi_r = [0, 0, 1]
    pi_v = [0.2, 0.8, 0]
    
    a = 600
    d = 1000
    
    data_point = DataPoint(a, d, pi_r, pi_v, mu)
    
    collapsed_sampler = CollapsedSampler()
    semi_collapsed_sampler = SemiCollapsedSampler()
    
    collapsed_results = collapsed_sampler.run(data_point, max_iters=100000, thin=100, burnin=50000)
    semi_collapsed_results = semi_collapsed_sampler.run(data_point, max_iters=1000000, thin=100, burnin=50000)
    
    from math import exp
    
    import matplotlib.pyplot as plot
    
    n = 1000
    likelihood = CollapsedLikelihood(data_point)
    
    x = []
    y = []
    
    norm_const = 0
    
    for i in range(1, n):
        x.append(i / n)
        y.append(exp(likelihood.compute_log_likelihood(x[-1])))
        
        norm_const += 1 / n * y[-1]
    
    y = [y_i / norm_const for y_i in y]
    
    plot.plot(x, y)
    plot.hist(collapsed_results['phi'], normed=True, bins=100, alpha=0.5)
    plot.hist(semi_collapsed_results['phi'], normed=True, bins=100, alpha=0.5)
    plot.show()

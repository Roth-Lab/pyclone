from __future__ import division

import bisect
import functools
import random

from collections import OrderedDict
from math import exp, log, lgamma as log_gamma, isinf

#=======================================================================================================================
# Log space normalisation.
#=======================================================================================================================
def log_sum_exp(log_X):
    '''
    Given a list of values in log space, log_X. Compute exp(log_X[0] + log_X[1] + ... log_X[n])
    
    Numerically safer than naive method.
    '''
    max_exp = max(log_X)
    
    if isinf(max_exp):
        return max_exp
    
    total = 0

    for x in log_X:
        total += exp(x - max_exp)
    
    return log(total) + max_exp

def log_space_normalise(log_X):
    normalised_log_X = []
    
    log_norm_const = log_sum_exp(log_X)
    
    for x in log_X:
        normalised_log_X.append(x - log_norm_const)
    
    return normalised_log_X

#=======================================================================================================================
# Distribution related code.
#=======================================================================================================================
def log_factorial(n):
    return log_gamma(n + 1)

def log_binomial_coefficient(n, k):
    return log_factorial(n) - log_factorial(k) - log_factorial(n - k)

def log_binomial_likelihood(x, n, mu):    
    if mu == 0:
        if x == 0:
            return 0
        else:
            return float('-inf')
    
    elif mu == 1:
        if x == n:
            return 0
        else:
            return float('-inf')
    
    else:    
        return x * log(mu) + (n - x) * log(1 - mu)

def log_binomial_pdf(x, n, mu):
    return log_binomial_coefficient(n, x) + x * log(mu) + (n - x) * log(1 - mu)

def log_beta(a, b):
    return log_gamma(a) + log_gamma(b) - log_gamma(a + b)

def log_beta_pdf(x, a, b):
    if x == 0 or x == 1:
        return float('-inf')    
    
    return -log_beta(a, b) + (a - 1) * log(x) + (b - 1) * log(1 - x)

#=======================================================================================================================
# Random variate generators
#=======================================================================================================================
def discrete_rvs(p):
    choices = range(len(p))
    
    cum_dist = [p[0]]
    
    for i in range(1, len(p)):
        cum_dist.append(cum_dist[i - 1] + p[i])
    
    x = random.random() * cum_dist[-1]
    
    return choices[bisect.bisect(cum_dist, x)]

def bernoulli_rvs(p):
    u = random.random()
    
    if u <= p:
        return 1
    else:
        return 0
    
#=======================================================================================================================
# Integration
#=======================================================================================================================
class Integrator(object):
    def __init__(self, a=0, b=1, mesh_size=100):
        self.a = a
        self.b = b
        self.mesh_size = mesh_size
        
        self.step_size = (b - a) / mesh_size
        
        self.knots = [i * self.step_size + a for i in range(0, mesh_size + 1)]
    
class SimpsonsRuleIntegrator(Integrator):
    def __init__(self, a=0, b=1, mesh_size=100):
        if mesh_size % 2 != 0:
            raise Exception("Mesh size for Simpson's rule must be an even number.")
        
        Integrator.__init__(self, a, b, mesh_size)

    def log_integrate(self, log_f):
        log_total = []
        
        # First and last terms.
        log_total.append(log_f(self.knots[0]))
        log_total.append(log_f(self.knots[-1]))
        
        four_total = []
        
        for i in range(1, self.mesh_size, 2):
            four_total.append(log_f(self.knots[i]))
        
        log_total.append(log(4) + log_sum_exp(four_total))
        
        two_total = []
        
        for i in range(2, self.mesh_size - 1, 2):
            two_total.append(log_f(self.knots[i]))
        
        log_total.append(log(2) + log_sum_exp(two_total))
  
        return log(self.step_size) - log(3) + log_sum_exp(log_total)

#=======================================================================================================================
# Function caching
#=======================================================================================================================
class memoized(object):
    def __init__(self, func, cache_size=10000):
        self.func = func
        
        self.cache = OrderedDict()

        self.cache_size = cache_size        
    
    def __call__(self, *args):
        if args in self.cache:
            value = self.cache[args]
        else:
            value = self.func(*args)
            
            self.cache[args] = value
        
        if len(self.cache) > self.cache_size:
            self.cache.popitem(last=False)
        
        return value

    def __get__(self, obj, objtype):
        '''Support instance methods.'''
        
        return functools.partial(self.__call__, obj)
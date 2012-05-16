from __future__ import division

import bisect
import random

from math import exp, log, lgamma as log_gamma

#=======================================================================================================================
# Log space normalisation.
#=======================================================================================================================
def log_sum_exp(log_X):
    '''
    Given a list of values in log space, log_X. Compute exp(log_X[0] + log_X[1] + ... log_X[n])
    
    Numerically safer than naive method.
    '''
    max_exp = log_X[0]
 
    for x in log_X:
        if max_exp < x:
            max_exp = x

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
    return x * log(mu) + (n - x) * log(1 - mu)

def log_binomial_pdf(x, n, mu):
    return log_binomial_coefficient(n, x) + x * log(mu) + (n - x) * log(1 - mu)

def log_beta(a, b):
    return log_gamma(a) + log_gamma(b) - log_gamma(a + b)

def log_beta_pdf(x, a, b):
    return log_beta(a, b) + (a - 1) * log(x) + (b - 1) * log(1 - x)

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

def histogram(values, num_bins=100, min_value=0, max_value=1, normalise=False):
    '''
    Simple histogram implementation. Intervals are evenly spaced half open and left inclusive.
    
    For example 0 would be in bin [0, 0.1) and 0.1 would be in bin [0.1,0.2).
    
    The right end-point is an exception so that if max_value=1 then the last bin is [0.9, 1].
    
    Args:
        values : (list) Values to create histogram from.
    
    Kwargs:
        num_bins : (int) Number of bins to use.
        min_value : (int) Left endpoint of histogram.
        max_value : (int) Right endpoint of histogram.
        normalise : (bool) If False then the integer counts of items in each bins is returned. 
                           If True returns the normalised histogram such that the integral sums to 1.
    '''     
    counts = [0 for _ in range(num_bins)]
    
    bin_width = (max_value - min_value) / num_bins
    
    num_endpoints = num_bins + 1
    
    endpoints = [(x * bin_width + min_value) for x in range(num_endpoints)]
    
    for x in values:
        bin_index = bisect.bisect_right(endpoints, x) - 1 
        
        # Corner case at right end-point
        if x == endpoints[-1]:
            counts[-1] += 1
        # Corner case when values outside endpoints. 
        elif x < endpoints[0] or x > endpoints[-1]:
            continue
        else:
            counts[bin_index] += 1
    
    if normalise == True:
        total = sum(counts)
        
        denom = bin_width * total
        
        hist = [x / denom for x in counts]
    else:
        hist = counts
    
    bin_centres = []
    
    for left, right in zip(endpoints[:-1], endpoints[1:]):
        bin_centres.append((left + right) / 2)
         
    
    return bin_centres, hist    

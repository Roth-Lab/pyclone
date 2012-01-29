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

def binomial_rvs(mu, d):
    '''
    Naive method based on representation of binomial as sum of Bernoulli variables.
    '''
    a = 0
    
    for _ in range(d):
        a += bernoulli_rvs(mu)
    
    return a

def poisson_rvs(expected_value):
    '''
    Inverse transform method from "Simulation" by S. Ross pg 51. Computation done in log space to avoid numerical
    issues.
    '''
    i = 0
    log_p = -expected_value 
    log_F = log_p
    
    log_lambda = log(expected_value)
    
    U = random.random()
    
    if U != 0:
        log_U = log(U)
    else:
        log_U = float('-inf')
    
    while log_F <= log_U:
        log_p = log_lambda + log_p - log(i + 1)
        log_F = log_sum_exp((log_F, log_p))
        
        i += 1
    
    return i    

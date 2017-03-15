'''
Created on 15 Mar 2017

@author: Andrew Roth
'''
from math import lgamma as log_gamma

import numpy as np

try:
    from numba import jit

except ImportError:
    def identity_decorator(*args, **kwargs):
        if len(args) == 1 and hasattr(args[0], '__call__'):
            return args[0]

        else:
            def _f(f):
                return f
            return _f

    jit = identity_decorator


@jit(cache=True, nopython=True)
def log_beta(a, b):
    if a <= 0 or b <= 0:
        return -np.inf

    return log_gamma(a) + log_gamma(b) - log_gamma(a + b)


@jit(cache=True, nopython=True)
def log_beta_binomial_likelihood(x, n, a, b):
    return log_beta(a + x, b + n - x) - log_beta(a, b)


@jit(cache=True, nopython=True)
def log_binomial_likelihood(x, n, p):
    if p == 0:
        if x == 0:
            return 0
        else:
            return -np.inf

    if p == 1:
        if x == n:
            return 0
        else:
            return -np.inf

    return x * np.log(p) + (n - x) * np.log(1 - p)


@jit(cache=True, nopython=True)
def log_sum_exp(log_X):
    '''
    Given a list of values in log space, log_X. Compute exp(log_X[0] + log_X[1] + ... log_X[n])

    Numerically safer than naive method.
    '''
    max_exp = np.max(log_X)

    if np.isinf(max_exp):
        return max_exp

    total = 0

    for x in log_X:
        total += np.exp(x - max_exp)

    return np.log(total) + max_exp

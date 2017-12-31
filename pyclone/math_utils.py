'''
Created on 15 Mar 2017

@author: Andrew Roth
'''
from math import lgamma as log_gamma

import numba
import numpy as np


@numba.jit(nopython=True)
def log_pyclone_beta_binomial_pdf(data, f, s):
    t = data.t

    C = len(data.cn)

    population_prior = np.zeros(3)
    population_prior[0] = (1 - t)
    population_prior[1] = t * (1 - f)
    population_prior[2] = t * f

    ll = np.ones(C, dtype=np.float64) * np.inf * -1

    for c in range(C):
        e_vaf = 0

        norm_const = 0

        for i in range(3):
            e_cn = population_prior[i] * data.cn[c, i]

            e_vaf += e_cn * data.mu[c, i]

            norm_const += e_cn

        e_vaf /= norm_const

        a = e_vaf * s

        b = s - a

        ll[c] = data.log_pi[c] + log_beta_binomial_pdf(data.a + data.b, data.b, a, b)

    return log_sum_exp(ll)


@numba.jit(nopython=True)
def log_pyclone_binomial_pdf(data, f):
    t = data.t

    C = len(data.cn)

    population_prior = np.zeros(3)
    population_prior[0] = (1 - t)
    population_prior[1] = t * (1 - f)
    population_prior[2] = t * f

    ll = np.ones(C, dtype=np.float64) * np.inf * -1

    for c in range(C):
        e_vaf = 0

        norm_const = 0

        for i in range(3):
            e_cn = population_prior[i] * data.cn[c, i]

            e_vaf += e_cn * data.mu[c, i]

            norm_const += e_cn

        e_vaf /= norm_const

        ll[c] = data.log_pi[c] + log_binomial_pdf(data.a + data.b, data.b, e_vaf)

    return log_sum_exp(ll)


@numba.jit(cache=True, nopython=True)
def log_beta(a, b):
    if a <= 0 or b <= 0:
        return -np.inf

    return log_gamma(a) + log_gamma(b) - log_gamma(a + b)


@numba.jit(cache=True, nopython=True)
def log_factorial(x):
    return log_gamma(x + 1)


@numba.jit(cache=True, nopython=True)
def log_binomial_coefficient(n, x):
    return log_factorial(n) - log_factorial(x) - log_factorial(n - x)


@numba.jit(cache=True, nopython=True)
def log_beta_binomial_likelihood(n, x, a, b):
    return log_beta(a + x, b + n - x) - log_beta(a, b)


@numba.jit(nopython=True)
def log_beta_binomial_pdf(n, x, a, b):
    return log_binomial_coefficient(n, x) + log_beta_binomial_likelihood(n, x, a, b)


@numba.jit(cache=True, nopython=True)
def log_binomial_likelihood(n, x, p):
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


@numba.jit(nopython=True)
def log_binomial_pdf(n, x, p):
    return log_binomial_coefficient(n, x) + log_binomial_likelihood(n, x, p)


@numba.jit(cache=True, nopython=True)
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


@numba.jit(cache=True, nopython=True)
def log_normalize(log_p):
    return log_p - log_sum_exp(log_p)

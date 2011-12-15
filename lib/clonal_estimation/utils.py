from math import exp, log, lgamma as log_gamma

def log_factorial(n):
    return log_gamma(n + 1)

def log_binomial_coefficient(n, k):
    return log_factorial(n) - log_factorial(k) - log_factorial(n - k)

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

def binomial_pdf(x, n, mu):
    return exp(log_binomial_coefficient(n, x)) * (mu ** x) * (1 - mu) ** (n - x)

def log_binomial_pdf(x, n, mu):
    return log_binomial_coefficient(n, x) + x * log(mu) + (n - x) * log(1 - mu)

def log_beta(a, b):
    return log_gamma(a) + log_gamma(b) - log_gamma(a + b)

def log_beta_pdf(x, a, b):
    return log_beta(a, b) + (a - 1) * log(x) + (b - 1) * log(1 - x)

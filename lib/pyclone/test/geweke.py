'''
Created on 2012-09-18

@author: Andrew Roth
'''
from __future__ import division

from math import sqrt

def compare_trace(trace_1, trace_2, g):
    g_1 = [g(x) for x in trace_1]
    g_2 = [g(x) for x in trace_2]
    
    mean_1 = mean(g_1)
    mean_2 = mean(g_2)
    
    sigma_1 = variance(g_1)
    sigma_2 = variance(g_2)

    M_1 = len(g_1)
    M_2 = len(g_2)
    
    diff = (mean_1 - mean_2) 
    
    print mean_1, mean_2, sigma_1, sigma_2
    
    return diff / sqrt(sigma_1 / M_1 + sigma_2 / M_2)

def mean(x):
    return sum(x) / len(x)

def variance(x):
    ss = sum([x_i ** 2 for x_i in x])
    
    return ss / len(x) - mean(x) ** 2

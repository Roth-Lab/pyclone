'''
Created on 2012-09-18

@author: Andrew Roth
'''
import numpy as np

def compare_trace(trace_1, trace_2, g):
    g_1 = g(trace_1)
    g_2 = g(trace_2)
    
    mean_1 = g_1.mean()
    mean_2 = g_2.mean()
    
    sigma_1 = g_1.var()
    sigma_2 = g_2.var()
    
    M_1 = len(g_1)
    M_2 = len(g_2)
    
    diff = (mean_1 - mean_2) 
    
    return diff / np.sqrt(sigma_1 / M_1 + sigma_2 / M_2)
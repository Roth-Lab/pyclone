'''
Utility functions for post-processing the results of a PyClone analysis.

Created on 2013-02-08

@author: Andrew Roth
'''
import pandas as pd

def load_cellular_frequencies_trace(file_name, burnin, thin):
    return _load_trace(file_name, burnin, thin, float)

def load_cluster_labels_trace(file_name, burnin, thin):
    return _load_trace(file_name, burnin, thin, int)

def _load_trace(trace_file, burnin, thin, cast_func):
    '''
        Args:
            trace_file : (str) Path to file to load.
            burnin : (int) Number of samples from the begining of MCMC chain to discard.
            thin : (int) Number of samples to skip when building trace.
            cast_func : (function) A function to cast data from string to appropriate type i.e. int, float
    '''        
    trace = pd.read_csv(trace_file, compression='bz2', sep='\t')
    
    trace = trace.iloc[burnin::thin]
    
    return trace.astype(cast_func)

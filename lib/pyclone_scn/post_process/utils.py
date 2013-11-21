'''
Utility functions for post-processing the results of a PyClone analysis.

Created on 2013-02-08

@author: Andrew Roth
'''
import bz2
import StringIO

try:
    import pandas as pd
except:
    raise Exception("The multi sample plotting module requires the pandas package. See http://http://pandas.pydata.org.")

def load_trace(trace_file, burnin, thin):
    '''
        Args:
            trace_file : (str) Path to file to load.
            burnin : (int) Number of samples from the begining of MCMC chain to discard.
            thin : (int) Number of samples to skip when building trace.
    '''
    fh = open(trace_file)
    
    data = StringIO.StringIO(bz2.decompress(fh.read()))
    
    data = pd.read_csv(data, sep='\t')
    
    data = data[burnin::thin]
    
    return data.T
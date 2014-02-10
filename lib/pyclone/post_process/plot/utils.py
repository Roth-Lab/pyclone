'''
Created on 2012-08-20

@author: Andrew Roth
'''
import csv
csv.field_size_limit(10000000)

def load_trace(trace_file, burnin, thin, cast_func):
    '''
        Args:
            trace_file : (str) Path to file to load.
            burnin : (int) Number of samples from the begining of MCMC chain to discard.
            thin : (int) Number of samples to skip when building trace.
            cast_func : (function) A function to cast data from string to appropriate type i.e. int, float
    '''    
    trace = {}
    
    reader = csv.DictReader(open(trace_file), delimiter='\t')
    
    for row in reader:
        gene = row['gene']
        
        trace[gene] = [cast_func(x) for x in row['trace'].split(',')][burnin::thin]
    
    return trace
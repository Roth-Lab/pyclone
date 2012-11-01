'''
Created on 2012-02-08

@author: Andrew Roth
'''
from math import lgamma as log_gamma

from pyclone.sampler import DirichletProcessSampler, PyCloneData
from pyclone.trace import TraceDB

import csv
import shutil

def run_dp_model(args):
    '''
    Run a fresh instance of the DP model.
    '''
    data, genes = load_pyclone_data(args.in_file)
    
    trace_db = TraceDB(args.out_dir, genes)
    
    try:
        sampler = DirichletProcessSampler(args.tumour_content, alpha=args.concentration)
    except:
        trace_db.close()
        
        shutil.rmtree(args.out_dir)
        
        raise
    
    sampler.sample(data, trace_db, num_iters=args.num_iters)

    trace_db.close()

def load_pyclone_data(file_name):
    '''
    Load data from PyClone formatted input file.
    '''
    data = []
    genes = []
    
    reader = csv.DictReader(open(file_name), delimiter='\t')

    for row in reader:
        genes.append(row['mutation'])
        
        a = int(row['a'])
        
        d = int(row['d'])
        
        mu_r = [float(x) for x in row['mu_r'].split(',')]
        mu_v = [float(x) for x in row['mu_v'].split(',')]
        
        delta_r = [float(x) for x in row['delta_r'].split(',')]
        delta_v = [float(x) for x in row['delta_v'].split(',')]
        
        log_pi_r = get_log_mix_weights(delta_r)
        log_pi_v = get_log_mix_weights(delta_v)
        
        data.append(PyCloneData(a, d, tuple(mu_r), tuple(mu_v), tuple(log_pi_r), tuple(log_pi_v)))

    return data, genes

def get_log_mix_weights(delta):
    log_denominator = log_gamma(sum(delta) + 1)
    
    log_mix_weights = []
    
    for i, d_i in enumerate(delta):
        log_numerator = log_gamma(d_i + 1)
        
        for j, d_j in enumerate(delta):
            if i != j:
                log_numerator += log_gamma(d_j)
        
        log_mix_weights.append(log_numerator - log_denominator)
    
    return log_mix_weights

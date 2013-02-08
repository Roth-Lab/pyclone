'''
Created on 2012-02-08

@author: Andrew Roth
'''
from __future__ import division

from math import log

import csv
import numpy as np
import shutil

from pyclone.sampler import DirichletProcessSampler, Data
from pyclone.trace import TraceDB

def run_dp_model(args):
    '''
    Run a fresh instance of the DP model.
    '''
    data = load_pyclone_data(args.in_file, args.error_rate, args.sampling_model)
    
    trace_db = TraceDB(args.out_dir, data.keys())
    
    try:
        sampler = DirichletProcessSampler(args.tumour_content,
                                          alpha=args.concentration,
                                          alpha_shape=args.concentration_prior_shape,
                                          alpha_rate=args.concentration_prior_rate)
    except:
        trace_db.close()
        
        shutil.rmtree(args.out_dir)
        
        raise
    
    sampler.sample(data.values(), trace_db, num_iters=args.num_iters)

    trace_db.close()

def load_pyclone_data(file_name, error_rate, sampling_model):
    '''
    Load data from PyClone formatted input file.
    '''
    data = {}
    
    reader = csv.DictReader(open(file_name), delimiter='\t')

    for row in reader:
        mutation = row['mutation']
        
        b = int(row['b'])
        
        d = int(row['d'])
        
        prior_weights = [float(x) for x in row['prior_weight'].split(',')]
            
        log_pi = get_log_mix_weights(prior_weights)
        
        mu_v = [float(x) for x in row['mu_v'].split(',')]        
        
        if sampling_model == 'fragment':
            cn_r = [float(x) for x in row['cn_r'].split(',')]
            
            cn_v = [float(x) for x in row['cn_v'].split(',')]
        
            data[mutation] = Data(b, d, error_rate, np.array(cn_r), np.array(cn_v), np.array(mu_v), np.array(log_pi))

    return data

def get_log_mix_weights(delta):
    pi = [x / sum(delta) for x in delta]
    
    return [log(x) for x in pi]

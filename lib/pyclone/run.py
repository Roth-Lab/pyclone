'''
Created on 2012-02-08

@author: Andrew Roth
'''
from __future__ import division

from math import log

import csv
import shutil

from pyclone.sampler import DirichletProcessSampler, DataPoint
from pyclone.trace import TraceDB

def run_dp_model(args):
    '''
    Run a fresh instance of the DP model.
    '''
    data = load_pyclone_data(args.in_file, args.error_rate)
    
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
    
    sampler.sample(data.values(), trace_db, num_iters=args.num_iters, seed=args.seed)

    trace_db.close()

def load_pyclone_data(file_name, error_rate):
    '''
    Load data from PyClone formatted input file.
    '''
    data = {}
    
    reader = csv.DictReader(open(file_name), delimiter='\t')

    for row in reader:
        mutation = row['mutation']
        
        b = int(row['b'])
        
        d = int(row['d'])
        
        weights = [float(x) for x in row['prior_weight'].split(',')]
        
        mu_v = [float(x) for x in row['mu_v'].split(',')]        
        
        cn_r = [float(x) for x in row['cn_r'].split(',')]
        
        cn_v = [float(x) for x in row['cn_v'].split(',')]
    
        data[mutation] = DataPoint(b, d, error_rate, cn_r, cn_v, mu_v, weights)

    return data
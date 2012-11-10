'''
Created on 2012-02-08

@author: Andrew Roth
'''
from __future__ import division

from math import log

import csv
import shutil

from pyclone.sampler import DirichletProcessSampler, PyCloneData
from pyclone.trace import TraceDB

def run_dp_model(args):
    '''
    Run a fresh instance of the DP model.
    '''
    data, mutations = load_pyclone_data(args.in_file)
    
    trace_db = TraceDB(args.out_dir, mutations)
    
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
    mutations = []
    
    reader = csv.DictReader(open(file_name), delimiter='\t')

    for row in reader:
        mutations.append(row['mutation'])
        
        a = int(row['a'])
        
        d = int(row['d'])
        
        mu_r = [float(x) for x in row['mu_r'].split(',')]
        mu_v = [float(x) for x in row['mu_v'].split(',')]
        
        delta_r = [float(x) for x in row['delta_r'].split(',')]
        delta_v = [float(x) for x in row['delta_v'].split(',')]
        
        log_pi_r = get_log_mix_weights(delta_r)
        log_pi_v = get_log_mix_weights(delta_v)
        
        data.append(PyCloneData(a, d, tuple(mu_r), tuple(mu_v), tuple(log_pi_r), tuple(log_pi_v)))

    return data, mutations

def get_log_mix_weights(delta):
    pi = [x / sum(delta) for x in delta]
    
    return [log(x) for x in pi]

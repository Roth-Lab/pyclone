'''
Created on 2012-02-08

@author: Andrew Roth
'''
from __future__ import division

from math import log

import csv
import shutil

from pyclone.densities import FragmentSampledData, CellSampledData, CellSampledDensity, FragmentSampledDensity
from pyclone.sampler import DirichletProcessSampler
from pyclone.trace import TraceDB

def run_dp_model(args):
    '''
    Run a fresh instance of the DP model.
    '''
    data, mutations = load_pyclone_data(args.in_file, args.sampling_model)
    
    trace_db = TraceDB(args.out_dir, mutations)
    
    if args.sampling_model == 'cell':
        cluster_density = CellSampledDensity()
    elif args.sampling_model == 'fragment':
        cluster_density = FragmentSampledDensity()
    
    try:
        sampler = DirichletProcessSampler(cluster_density,
                                          args.tumour_content,
                                          alpha=args.concentration,
                                          alpha_shape=args.concentration_prior_shape,
                                          alpha_rate=args.concentration_prior_rate)
    except:
        trace_db.close()
        
        shutil.rmtree(args.out_dir)
        
        raise
    
    sampler.sample(data, trace_db, num_iters=args.num_iters)

    trace_db.close()

def load_pyclone_data(file_name, sampling_model):
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
        
        if sampling_model == 'dna':
            cn_r = [int(x) for x in row['cn_r'].split(',')]
            cn_v = [int(x) for x in row['cn_v'].split(',')]
            
            data_point = FragmentSampledData(a,
                                             d,
                                             tuple(mu_r),
                                             tuple(mu_v),
                                             tuple(log_pi_r),
                                             tuple(log_pi_v),
                                             tuple(cn_r),
                                             tuple(cn_v))
        
        else:
            data_point = CellSampledData(a, d, tuple(mu_r), tuple(mu_v), tuple(log_pi_r), tuple(log_pi_v))
            
        data.append(data_point)

    return data, mutations

def get_log_mix_weights(delta):
    pi = [x / sum(delta) for x in delta]
    
    return [log(x) for x in pi]

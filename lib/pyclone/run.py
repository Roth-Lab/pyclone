'''
Created on 2012-02-08

@author: Andrew Roth
'''
from __future__ import division

import csv
import os
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

def cluster_trace(args):
    from pyclone.post_process.cluster import cluster_pyclone_trace
    
    pyclone_file = os.path.join(args.trace_dir, 'labels.tsv.bz2')
    
    print '''Clustering PyClone trace file {in_file} using {method} with a burnin of {burnin} and using every {thin}th sample'''.format(in_file=pyclone_file, method=args.method, burnin=args.burnin, thin=args.thin)    
    
    cluster_pyclone_trace(pyclone_file, args.out_file, args.method, args.burnin, args.thin)

def plot_cellular_frequencies(args):
    import pyclone.post_process.plot as plot
    
    pyclone_file = os.path.join(args.trace_dir, 'cellular_frequencies.tsv.bz2')
    
    print '''Plotting cellular frequencies from the PyClone trace file {in_file} with a burnin of {burnin} and using every {thin}th sample'''.format(in_file=pyclone_file, burnin=args.burnin, thin=args.thin)   
    
    plot.plot_cellular_frequencies(pyclone_file, args.out_file, args.burnin, args.thin) 
    
def plot_similarity_matrix(args):
    import pyclone.post_process.plot as plot
    
    pyclone_file = os.path.join(args.trace_dir, 'labels.tsv.bz2')
    
    print '''Plotting similarity matrix from the PyClone trace file {in_file} with a burnin of {burnin} and using every {thin}th sample'''.format(in_file=pyclone_file, burnin=args.burnin, thin=args.thin)   
    
    plot.plot_similarity_matrix(pyclone_file, args.out_file, args.burnin, args.thin) 

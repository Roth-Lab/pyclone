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
    

def build_prior_file(args):
    reader = csv.DictReader(open(args.in_file), delimiter='\t')
    
    extra_fields = [x for x in reader.fieldnames if x not in ['mutation', 'b', 'd']]
    
    writer_fields = ['mutation', 'b', 'd', 'cn_r', 'cn_v', 'mu_v', 'prior_weight'] + extra_fields
    
    writer = csv.DictWriter(open(args.out_file, 'w'), writer_fields, delimiter='\t')
    
    writer.writeheader()
    
    for row in reader:
        out_row = _get_prior_row(row, args.cn_r, args.mu_v, args.error_rate)
                
        for field_name in extra_fields:
            out_row[field_name] = row[field_name]
        
        writer.writerow(out_row)
    
def _get_prior_row(row, cn_r_method, mu_v_method, eps, cn_n=2):
    cn = int(row['cn'])
    
    mu_v = []
    
    if mu_v_method == 'single':
        mu_v.append(1 / cn)
    
    elif mu_v_method == 'all':
        mu_v.append(1 - eps)
    
    elif mu_v_method == 'vague':
        for num_var_alleles in range(1, cn):
            mu_v.append(num_var_alleles / cn)
        
        mu_v.append(1 - eps)
    
    if cn_r_method == 'normal':
        cn_r = [cn_n for _ in mu_v]
    
    elif cn_r_method == "variant":
        cn_r = [cn for _ in mu_v]
    
    elif cn_r_method == "vague":
        if cn == cn_n:
            cn_r = [cn for _ in mu_v]
        else:            
            cn_r = [cn_n for _ in mu_v] + [cn for _ in mu_v]
        
            mu_v = mu_v + mu_v
        
    cn_v = [cn for _ in mu_v]
    
    prior_weight = [1 for _ in mu_v]
    
    out_row = {
               'mutation' : row['mutation'],
               'b' : row['b'],
               'd' : row['d'],
               'cn_r' : list_to_csv(cn_r),
               'cn_v' : list_to_csv(cn_v),
               'mu_v' : list_to_csv(mu_v),
               'prior_weight' : list_to_csv(prior_weight)
               }
    
    return out_row

def list_to_csv(l):
    return ",".join([str(x) for x in l])

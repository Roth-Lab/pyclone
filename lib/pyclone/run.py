'''
Created on 2012-02-08

@author: Andrew Roth
'''
from __future__ import division

import csv
import random
import os
import yaml

from pyclone.config import get_mutation
from pyclone.ibmm import run_ibmm_analysis
from pyclone.igmm import run_igmm_analysis
from pyclone.utils import make_parent_directory

#=======================================================================================================================
# PyClone analysis
#=======================================================================================================================
def run_analysis(args):
    if args.seed is not None:
        random.seed(args.seed)
    
    fh = open(args.config_file)
    
    config = yaml.load(fh)
    
    fh.close()
    
    trace_dir = os.path.join(config['working_dir'], config['trace_dir'])
    
    alpha = config['concentration']['value']
    
    if 'prior' in config['concentration']:
        alpha_priors = config['concentration']['prior']
    else:
        alpha_priors = None
    
    num_iters = config['num_iters']
    
    density = config['density']
    
    if density == 'binomial':
        run_ibmm_analysis(args.config_file, trace_dir, num_iters, alpha, alpha_priors)
    
    elif density =='gaussian':
        run_igmm_analysis(args.config_file, trace_dir, num_iters, alpha, alpha_priors)
    

# def run_dp_model(args):
#     '''
#     Run a fresh instance of the DP model.
#     '''
#     data = load_pyclone_data(args.in_file, args.tumour_content)
#     
#     trace = DiskTrace(args.out_dir, 
#                       ['alpha', 'labels', 'x'], 
#                       column_names=data.keys(), 
#                       file_name_map={'x' : 'cellular_frequencies'})
#     
#     trace.open('w')
#     
#     try:
#         sampler = PyCloneSampler(alpha=args.concentration,
#                                  alpha_shape=args.concentration_prior_shape,
#                                  alpha_rate=args.concentration_prior_rate)
#     except:
#         trace.close()
#         
#         shutil.rmtree(args.out_dir)
#         
#         raise
#     
#     sampler.sample(data.values(), trace, num_iters=args.num_iters, seed=args.seed)
# 
#     trace.close()
# 
# def load_pyclone_data(file_name, tumour_content):
#     '''
#     Load data from PyClone formatted input file.
#     '''
#     data = OrderedDict()
#     
#     fh = open(file_name)
#     
#     config = yaml.load(fh)
#     
#     fh.close()
#     
#     error_rate = config['error_rate']
# 
#     for mutation_dict in config['mutations']:
#         mutation = load_mutation_from_dict(mutation_dict)
# 
#         data[mutation.id] = PyCloneData(mutation.ref_counts, 
#                                         mutation.var_counts, 
#                                         mutation.states, 
#                                         tumour_content, 
#                                         error_rate)
#     
#     return data
# 
# #=======================================================================================================================
# # DP Analysis code
# #=======================================================================================================================
# def run_dirichlet_process_analysis(args):
#     '''
#     Runs a genotype naive analysis using an infinte Gaussian, Binomial, or Beta-Binomial mixture model.
#     '''
#     if args.density in ['binomial', 'beta_binomial']:
#         data = load_count_data(args.in_file)
#=======================================================================================================================
# Input file code
#=======================================================================================================================
def build_mutation_file(args):
    config = {}
    
    reader = csv.DictReader(open(args.in_file), delimiter='\t')
    
    config['mutations'] = []

    for row in reader:
        mutation_id = row['mutation_id']
        
        ref_counts = int(row['ref_counts'])
        
        var_counts = int(row['var_counts'])
        
        normal_cn = int(row['normal_cn'])
        
        minor_cn = int(row['minor_cn'])
        
        major_cn = int(row['major_cn'])

        mutation = get_mutation(mutation_id, 
                                ref_counts, 
                                var_counts, 
                                normal_cn, 
                                minor_cn, 
                                major_cn, 
                                args.ref_prior, 
                                args.var_prior)

        config['mutations'].append(mutation.to_dict())
    
    make_parent_directory(args.out_file)
    
    fh = open(args.out_file, 'w')
    
    yaml.dump(config, fh)
    
    fh.close()

#=======================================================================================================================
# Post processing code
#=======================================================================================================================
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
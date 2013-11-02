'''
Created on 2012-02-08

@author: Andrew Roth
'''
from __future__ import division

import csv
import os
import random
import yaml

from pyclone.config import get_mutation
from pyclone.model import run_pyclone_analysis

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
    
    alpha = args.alpha_init
    
    if 'prior' in config['concentration']:
        alpha_priors = config['concentration']['prior']
    else:
        alpha_priors = None
    
    num_iters = config['num_iters']

    run_pyclone_analysis(args.config_file, trace_dir, num_iters, alpha, alpha_priors)

#=======================================================================================================================
# Input file code
#=======================================================================================================================
def get_mutations(file_name):
    mutations = []
    
    reader = csv.DictReader(open(file_name), delimiter='\t')

    for row in reader:
        mutation_id = row['mutation_id']
        
        ref_counts = int(row['ref_counts'])
        
        var_counts = int(row['var_counts'])
        
        normal_cn = int(row['normal_cn'])
        
        minor_cn = int(row['minor_cn'])
        
        major_cn = int(row['major_cn'])
        
        cn_prevalence = float(row['cn_prevalence'])

        mutation = get_mutation(mutation_id,
                                ref_counts,
                                var_counts,
                                normal_cn,
                                minor_cn,
                                major_cn,
                                cn_prevalence)

        mutations.append(mutation)

#=======================================================================================================================
# Post processing code
#=======================================================================================================================
def build_multi_sample_table(args):
    from pyclone.post_process.plot.multi_sample import load_multi_sample_table
    
    table = load_multi_sample_table(args.config_file, args.prevalence, args.clustering_method, args.burnin, args.thin)

    table.to_csv(args.out_file, sep='\t')
        
def cluster_trace(args):
    from pyclone.post_process.cluster import write_pyclone_cluster_file
    
    config = _load_yaml_config(args.config_file)
    
    labels_file = os.path.join(config['working_dir'], config['trace_dir'], 'labels.tsv.bz2')
    
    print '''Clustering PyClone trace file {in_file} using {method} with a burnin of {burnin} and using every {thin}th sample'''.format(in_file=labels_file,
                                                                                                                                        method=args.method,
                                                                                                                                        burnin=args.burnin,
                                                                                                                                        thin=args.thin)    
    
    write_pyclone_cluster_file(labels_file, args.out_file, args.method, args.burnin, args.thin)

def plot_cellular_frequencies(args):
    import pyclone.post_process.plot as plot
    
    config = _load_yaml_config(args.config_file)
    
    trace_dir = os.path.join(config['working_dir'], config['trace_dir'])
    
    for sample_id in config['samples']:
        file_name = os.path.join(trace_dir, '{0}.cellular_frequencies.tsv.bz2'.format(sample_id))
        
        print '''Plotting cellular frequencies from the PyClone trace file {in_file} with a burnin of {burnin} and using every {thin}th sample'''.format(in_file=file_name,
                                                                                                                                                         burnin=args.burnin,
                                                                                                                                                         thin=args.thin)   
        
        out_file = os.path.basename(file_name).replace('tsv.bz2', 'pdf')
        
        out_file = os.path.join(args.out_dir, out_file)
        
        plot.plot_cellular_frequencies(file_name, out_file, args.burnin, args.thin) 
    
def plot_similarity_matrix(args):
    import pyclone.post_process.plot as plot
    
    config = _load_yaml_config(args.config_file)
    
    labels_file = os.path.join(config['working_dir'], config['trace_dir'], 'labels.tsv.bz2')
    
    print '''Plotting similarity matrix from the PyClone trace file {in_file} with a burnin of {burnin} and using every {thin}th sample'''.format(in_file=labels_file,
                                                                                                                                                  burnin=args.burnin,
                                                                                                                                                  thin=args.thin)   
    
    plot.plot_similarity_matrix(labels_file, args.out_file, args.burnin, args.thin)
    
def plot_multi_sample(args):
    from pyclone.post_process.plot.multi_sample import plot_clusters, plot_mutations
    
    if args.separate_lines:
        plot_mutations(args.config_file,
                       args.plot_file,
                       args.prevalence,
                       args.clustering_method,
                       args.burnin,
                       args.thin)
        
    else:
        plot_clusters(args.config_file,
                      args.plot_file,
                      args.prevalence,
                      args.clustering_method,
                      args.burnin,
                      args.thin)

def _load_yaml_config(file_name):
    fh = open(file_name)
    
    config = yaml.load(fh)
    
    fh.close()
    
    return config        

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
from pyclone.pyclone_beta_binomial import run_pyclone_beta_binomial_analysis
from pyclone.pyclone_binomial import run_pyclone_binomial_analysis
from pyclone.utils import make_parent_directory

import pyclone.post_process as post_process
import pyclone.post_process.plot as plot

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
                
    if density == 'pyclone_beta_binomial':
        run_pyclone_beta_binomial_analysis(args.config_file, trace_dir, num_iters, alpha, alpha_priors)   
    
    elif density == 'pyclone_binomial':
        run_pyclone_binomial_analysis(args.config_file, trace_dir, num_iters, alpha, alpha_priors)
    
    else:
        raise Exception('{0} is not a valid density for PyClone.'.format(density))   

#=======================================================================================================================
# Input file code
#=======================================================================================================================
def build_mutations_file(args):
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
def write_multi_sample_table(args):
    table = post_process.load_multi_sample_table(args.config_file, args.burnin, args.thin, old_style=args.old_style)
    
    table.to_csv(args.out_file, index=False, sep='\t')
    
def write_clusters_trace(args):
    with open(args.config_file) as fh:
        config = yaml.load(fh)
        
    labels_file = os.path.join(config['working_dir'], config['trace_dir'], 'labels.tsv.bz2')
    
    print '''Clustering PyClone trace file {in_file} with a burnin of {burnin} and using every {thin}th sample'''.format(in_file=labels_file,
                                                                                                                         burnin=args.burnin, 
                                                                                                                         thin=args.thin)    
    labels = post_process.cluster_pyclone_trace(labels_file, args.burnin, args.thin)
    
    labels.to_csv(args.out_file, index=False, sep='\t')
   
def plot_cellular_prevalence_posteriors(args):
    with open(args.config_file) as fh:
        config = yaml.load(fh)
        
    trace_dir = os.path.join(config['working_dir'], config['trace_dir'])
    
    labels_trace_file = os.path.join(trace_dir, 'labels.tsv.bz2')
    
    for sample_id in config['samples']:    
        prevalence_trace_file = os.path.join(trace_dir, '{0}.cellular_frequencies.tsv.bz2'.format(sample_id))
        
        print '''Plotting cellular frequencies from the PyClone trace file {in_file} with a burnin of {burnin} and using every {thin}th sample'''.format(in_file=prevalence_trace_file, 
                                                                                                                                                         burnin=args.burnin, 
                                                                                                                                                         thin=args.thin)   
        
        out_file = os.path.basename(prevalence_trace_file).replace('tsv.bz2', args.file_format)
        
        out_file = os.path.join(args.out_dir, out_file)
        
        plot.plot_cellular_prevalence_posteriors(labels_trace_file, prevalence_trace_file, out_file, args.burnin, args.thin) 
    
def plot_similarity_matrix(args):
    with open(args.config_file) as fh:
        config = yaml.load(fh)
    
    labels_file = os.path.join(config['working_dir'], config['trace_dir'], 'labels.tsv.bz2')
    
    print '''Plotting similarity matrix from the PyClone trace file {in_file} with a burnin of {burnin} and using every {thin}th sample'''.format(in_file=labels_file, 
                                                                                                                                                  burnin=args.burnin, 
                                                                                                                                                  thin=args.thin)   
    
    plot.plot_similarity_matrix(labels_file, args.out_file, args.burnin, args.thin)
    
def plot_multi_sample(args):
    plot.plot_multi_sample_parallel_coordinates(
        args.config_file, 
        args.plot_file, 
        args.y_value,  
        burnin=args.burnin, 
        thin=args.thin,
        samples=args.samples,
        separate_lines=args.separate_lines
    )

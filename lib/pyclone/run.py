'''
Created on 2012-02-08

@author: Andrew Roth
'''
from __future__ import division

import csv
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
    
    alpha = config['concentration']['value']
    
    if 'prior' in config['concentration']:
        alpha_priors = config['concentration']['prior']
    else:
        alpha_priors = None
    
    num_iters = config['num_iters']
    
    density = config['density']
                
    if density == 'pyclone_beta_binomial':
        run_pyclone_beta_binomial_analysis(
            args.config_file, 
            num_iters, 
            alpha, 
            alpha_priors
        )   
    
    elif density == 'pyclone_binomial':
        run_pyclone_binomial_analysis(
            args.config_file, 
            num_iters, 
            alpha, 
            alpha_priors
        )
    
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
    labels = post_process.cluster_pyclone_trace(args.config_file, args.burnin, args.thin)
    
    labels.to_csv(args.out_file, index=False, sep='\t')
   
def plot_cellular_prevalence_posteriors(args):
    plot.plot_cellular_prevalence_posteriors(
        args.config_file, 
        args.file_format,
        args.out_dir, 
        args.burnin, 
        args.thin
    ) 

def plot_similarity_matrix(args):
    plot.plot_similarity_matrix(
        args.config_file, 
        args.out_file, 
        args.burnin, 
        args.thin
    )
    
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

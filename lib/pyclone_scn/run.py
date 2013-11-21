'''
Created on 2012-02-08

@author: Andrew Roth
'''
from __future__ import division

import os
import random
import yaml

from pyclone_scn.model import run_pyclone_analysis
from pyclone_scn.utils import make_directory

#=======================================================================================================================
# PyClone analysis
#=======================================================================================================================
def run_analysis(args):
    if args.seed is not None:
        random.seed(args.seed)
    
    config = _load_yaml_config(args.config_file)
    
    make_directory(config['out_dir'])
    
    trace_dir = os.path.join(config['out_dir'], 'trace')
    
    alpha = args.alpha_init
    
    alpha_priors = {'shape' : args.alpha_shape, 'rate' : args.alpha_rate}

    num_iters = args.num_iters

    run_pyclone_analysis(args, args.config_file, trace_dir, num_iters, alpha, alpha_priors)

#=======================================================================================================================
# Post processing code
#=======================================================================================================================
def write_results_table(args):
    from pyclone_scn.post_process.plot.multi_sample import load_results_table
    
    table = load_results_table(args.config_file, args.burnin, args.thin)

    table.to_csv(args.out_file, index_label='mutation_id', float_format='%.4f', sep='\t')
 
def plot_cellular_prevalences(args):
    import pyclone_scn.post_process.plot as plot
    
    config = _load_yaml_config(args.config_file)
    
    trace_dir = os.path.join(config['out_dir'], 'trace')
    
    make_directory(args.out_dir)
    
    for entry in config['samples']:
        sample_id = entry['sample_id']
        
        file_name = os.path.join(trace_dir, '{0}.cellular_prevalences.tsv.bz2'.format(sample_id))
        
        print '''Plotting cellular frequencies from the PyClone trace file {in_file} with a burnin of {burnin} and using every {thin}th sample'''.format(in_file=file_name,
                                                                                                                                                         burnin=args.burnin,
                                                                                                                                                         thin=args.thin)   
        
        out_file = os.path.basename(file_name).replace('tsv.bz2', 'pdf')
        
        out_file = os.path.join(args.out_dir, out_file)
        
        plot.plot_cellular_frequencies(file_name, out_file, args.burnin, args.thin) 
    
def plot_similarity_matrix(args):
    import pyclone_scn.post_process.plot as plot
    
    config = _load_yaml_config(args.config_file)
    
    labels_file = os.path.join(config['out_dir'], 'trace', 'labels.tsv.bz2')
    
    print '''Plotting similarity matrix from the PyClone trace file {in_file} with a burnin of {burnin} and using every {thin}th sample'''.format(in_file=labels_file,
                                                                                                                                                  burnin=args.burnin,
                                                                                                                                                  thin=args.thin)   
    
    plot.plot_similarity_matrix(labels_file, args.out_file, args.burnin, args.thin)
    
def plot_multi_sample(args):
    from pyclone_scn.post_process.plot.multi_sample import plot_clusters, plot_mutations
    
    if args.separate_lines:
        plot_mutations(args.config_file,
                       args.plot_file,
                       args.prevalence,
                       args.burnin,
                       args.thin)
        
    else:
        plot_clusters(args.config_file,
                      args.plot_file,
                      args.prevalence,
                      args.burnin,
                      args.thin)

def _load_yaml_config(file_name):
    fh = open(file_name)
    
    config = yaml.load(fh)
    
    fh.close()
    
    return config        

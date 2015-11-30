'''
Created on 2012-02-08

@author: Andrew Roth
'''
from __future__ import division

from collections import OrderedDict

import csv
import os
import random
import yaml

from pyclone.config import get_mutation
from pyclone.pyclone_beta_binomial import run_pyclone_beta_binomial_analysis
from pyclone.pyclone_binomial import run_pyclone_binomial_analysis
from pyclone.utils import make_directory, make_parent_directory

import pyclone.paths as paths
import pyclone.post_process as post_process
import pyclone.post_process.plot as plot

#=======================================================================================================================
# PyClone analysis
#=======================================================================================================================
def run_analysis_pipeline(args):
    make_directory(args.working_dir)
    
    make_directory(os.path.join(args.working_dir, 'yaml'))
    
    mutations_files = OrderedDict()
    
    tumour_contents = {}
    
    for i, in_file in enumerate(args.in_files):
        if args.samples is not None:
            sample_id = args.samples[i]
            
        else:
            sample_id = os.path.splitext(os.path.basename(in_file))[0]
        
        mutations_files[sample_id] = os.path.join(args.working_dir, 'yaml', '{0}.yaml'.format(sample_id))
        
        _build_mutations_file(
            in_file, 
            mutations_files[sample_id], 
            args.var_prior
        )
        
        if args.tumour_contents is not None:
            tumour_contents[sample_id] = args.tumour_contents[i]
        
        else:
            tumour_contents[sample_id] = 1.0
    
    config_file = os.path.join(args.working_dir, 'config.yaml')
    
    _write_config_file(
        config_file,
        args.density, 
        mutations_files, 
        args.num_iters, 
        tumour_contents, 
        args.working_dir
    )
    
    _run_analysis(config_file, args.seed)
    
    table_file = os.path.join(args.working_dir, 'results.tsv')
    
    _write_multi_sample_table(
        config_file, 
        table_file, 
        args.burnin, 
        args.thin, 
        False
    )
    
    clusters_file = os.path.join(args.working_dir, 'clusters.tsv')
    
    _write_clusters_trace(
        config_file, 
        clusters_file, 
        args.burnin, 
        args.thin
    )
    
    cluster_posteriors_file = os.path.join(args.working_dir, 'cluster_posteriors.tsv')
    
    _write_cluster_posteriors_table(
        config_file, 
        cluster_posteriors_file, 
        args.burnin, 
        args.thin, 
        101)
    
    plots_dir = os.path.join(args.working_dir, 'plots')
    
    cellular_prevalence_posteriors_dir = os.path.join(plots_dir, 'cellular_prevalence_posteriors')
    
    make_directory(cellular_prevalence_posteriors_dir)
    
    _plot_cellular_prevalence_posteriors(
        config_file, 
        cellular_prevalence_posteriors_dir, 
        args.burnin, 
        args.thin, 
        args.plot_file_format)
    
    sim_mat_file = os.path.join(plots_dir, 'similarity_matrix.{0}'.format(args.plot_file_format))
    
    _plot_similarity_matrix(
        config_file, 
        sim_mat_file, 
        args.burnin, 
        args.thin
    )
    
    cellular_prevalence_parallel_coordinates_file = os.path.join(
        plots_dir, 
        'cellular_prevalence_parallel_coordinates.{0}'.format(args.plot_file_format)
    )
    
    _plot_multi_sample(
        config_file, 
        cellular_prevalence_parallel_coordinates_file, 
        args.burnin, 
        args.thin, 
        mutations_files.keys(), 
        False, 
        'cellular_prevalence'
    )
    
    cluster_posteriors_plot_file = os.path.join(
        plots_dir,
        'cluster_cellular_prevalence_posteriors.{0}'.format(args.plot_file_format)
    )
    
    _plot_cluster_posteriors(
        config_file, 
        cluster_posteriors_plot_file, 
        args.burnin, 
        args.thin, 
        101, 
        mutations_files.keys()
    )
    
    vaf_parallel_coordinates_file = os.path.join(
        plots_dir, 
        'vaf_parallel_coordinates.{0}'.format(args.plot_file_format)
    )
    
    _plot_multi_sample(
        config_file, 
        vaf_parallel_coordinates_file, 
        args.burnin, 
        args.thin, 
        mutations_files.keys(), 
        True, 
        'variant_allele_frequency'
    )

def _write_config_file(config_file, density, mutations_files, num_iters, tumour_contents, working_dir):
    config = {}
    
    config['num_iters'] = num_iters
        
    config['base_measure_params'] = {'alpha' : 1, 'beta' : 1}
    
    config['concentration'] = {
        'value' : 1.0, 
        'prior' :{
            'shape' : 1.0,
            'rate' : 0.001
            }
    }
    
    config['density'] = density
    
    if density == 'pyclone_beta_binomial':
        config['beta_binomial_precision_params'] = {
            'value' : 1000,
            'prior' : {
                'shape' : 1.0,
                'rate' : 0.001
            },
            'proposal' : {'precision' : 0.01}
        }
    
    config['working_dir'] = os.path.abspath(working_dir)
    
    config['trace_dir'] = 'trace'
    
    config['samples'] = {}
    
    for sample_id in mutations_files:
        config['samples'][sample_id] = {
            'mutations_file' : mutations_files[sample_id],
            'tumour_content' : {
                'value' : tumour_contents[sample_id] 
            },
            'error_rate' : 0.001
        }
    
    with open(config_file, 'w') as fh:
        yaml.dump(config, fh, default_flow_style=False)

def run_analysis(args):
    _run_analysis(args.config_file, args.seed)

def _run_analysis(config_file, seed):
    if seed is not None:
        random.seed(seed)
    
    config = paths.load_config(config_file)
    
    alpha = config['concentration']['value']
    
    if 'prior' in config['concentration']:
        alpha_priors = config['concentration']['prior']
    else:
        alpha_priors = None
    
    num_iters = config['num_iters']
    
    density = config['density']
                
    if density == 'pyclone_beta_binomial':
        run_pyclone_beta_binomial_analysis(
            config_file, 
            num_iters, 
            alpha, 
            alpha_priors
        )   
    
    elif density == 'pyclone_binomial':
        run_pyclone_binomial_analysis(
            config_file, 
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
    _build_mutations_file(
        args.in_file, 
        args.out_file,  
        args.var_prior
    )

def _build_mutations_file(in_file, out_file, var_prior):
    config = {}
    
    reader = csv.DictReader(open(in_file), delimiter='\t')
    
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
                                var_prior)

        config['mutations'].append(mutation.to_dict())
    
    make_parent_directory(out_file)
    
    fh = open(out_file, 'w')
    
    yaml.dump(config, fh)
    
    fh.close()

#=======================================================================================================================
# Post processing code
#=======================================================================================================================
def write_multi_sample_table(args):
    _write_multi_sample_table(
        args.config_file, 
        args.out_file, 
        args.burnin, 
        args.thin, 
        args.old_style
    )

def _write_multi_sample_table(config_file, out_file, burnin, thin, old_style):
    table = post_process.load_multi_sample_table(config_file, burnin, thin, old_style=old_style)
    
    table.to_csv(out_file, index=False, sep='\t')
    
def write_clusters_trace(args):
    _write_clusters_trace(
        args.config_file,
        args.out_file,
        args.burnin,
        args.thin)

def _write_clusters_trace(config_file, out_file, burnin, thin):
    labels = post_process.cluster_pyclone_trace(config_file, burnin, thin)
    
    labels.to_csv(out_file, index=False, sep='\t')

def write_cluster_posteriors_table(args):
    _write_cluster_posteriors_table(
        args.config_file, 
        args.out_file, 
        args.burnin, 
        args.thin, 
        args.mesh_size
    )

def _write_cluster_posteriors_table(config_file, out_file, burnin, thin, mesh_size):
    df = post_process.load_cluster_posteriors_table(
        config_file, 
        burnin=burnin, 
        thin=thin, 
        mesh_size=mesh_size
    )
    
    df.to_csv(out_file, index=False, sep='\t')

def plot_cellular_prevalence_posteriors(args):
    _plot_cellular_prevalence_posteriors(
        args.config_file, 
        args.out_dir, 
        args.burnin, 
        args.thin, 
        args.file_format)
    
def _plot_cellular_prevalence_posteriors(config_file, out_dir, burnin, thin, file_format):
    plot.plot_cellular_prevalence_posteriors(
        config_file, 
        file_format,
        out_dir, 
        burnin, 
        thin
    )
    
def plot_cluster_posteriors(args):
    _plot_cluster_posteriors(
        args.config_file, 
        args.plot_file, 
        args.burnin, 
        args.mesh_size,
        args.min_cluster_size,
        args.plot_type, 
        args.samples,
        args.thin, 
    )

def _plot_cluster_posteriors(config_file, plot_file, burnin, mesh_size, min_cluster_size, plot_type, samples, thin):
    if plot_type == 'density':
        plot.cluster_posteriors.density_plot(
            config_file, 
            plot_file, 
            burnin=burnin, 
            thin=thin, 
            mesh_size=mesh_size,
            min_cluster_size=min_cluster_size,
            samples=samples, 
        )
    
    elif plot_type == 'parallel_coordinates':
        plot.cluster_posteriors.parallel_coordinates(
            config_file, 
            plot_file, 
            burnin=burnin, 
            mesh_size=mesh_size,
            min_cluster_size=min_cluster_size,
            samples=samples, 
            thin=thin
        )

def plot_similarity_matrix(args):
    _plot_similarity_matrix(
        args.config_file, 
        args.plot_file, 
        args.burnin, 
        args.thin
    )

def _plot_similarity_matrix(config_file, plot_file, burnin, thin):
    plot.plot_similarity_matrix(
        config_file, 
        plot_file, 
        burnin, 
        thin
    )
    
def plot_multi_sample(args):
    _plot_multi_sample(
        args.config_file, 
        args.plot_file, 
        args.burnin, 
        args.thin, 
        args.samples, 
        args.separate_lines, 
        args.y_value
    )

def _plot_multi_sample(config_file, plot_file, burnin, thin, samples, separate_lines, y_value):
    plot.plot_multi_sample_parallel_coordinates(
        config_file,
        plot_file,
        y_value,
        burnin=burnin,
        thin=thin,
        samples=samples,
        separate_lines=separate_lines
    )


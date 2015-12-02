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
            args.prior
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
    
    tables_dir = os.path.join(args.working_dir, 'tables')
    
    make_directory(tables_dir)
    
    for table_type in ['cluster', 'loci']:
        out_file = os.path.join(tables_dir, '{0}.tsv'.format(table_type))
        
        _build_table(
            config_file, 
            out_file, 
            args.burnin, 
            args.mesh_size, 
            table_type, 
            args.thin
        )
        
    plots_dir = os.path.join(args.working_dir, 'plots')

    plots = [
        ('cluster', 'density'),
        ('cluster', 'parallel_coordinates'),
        ('cluster', 'scatter'),
        ('loci', 'density'),
        ('loci', 'parallel_coordinates'),
        ('loci', 'scatter'),
        ('loci', 'similarity_matrix'),
        ('loci', 'vaf_parallel_coordinates'),
        ('loci', 'vaf_scatter')
    ]
    
    for category, plot_type in plots:
        
        plot_file = os.path.join(plots_dir, category, '{0}.{1}'.format(plot_type, args.plot_file_format))
        
        make_parent_directory(plot_file)
        
        if category == 'cluster':
            
            _cluster_plot(
                config_file, 
                plot_file, 
                args.burnin, 
                args.mesh_size, 
                args.min_cluster_size, 
                plot_type, 
                args.samples, 
                args.thin
            )
            
        elif category == 'loci':
            
            _loci_plot(
                config_file, 
                plot_file, 
                plot_type, 
                args.burnin, 
                args.min_cluster_size, 
                args.samples, 
                args.thin
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
        args.prior
    )

def _build_mutations_file(in_file, out_file, prior):
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
                                prior)

        config['mutations'].append(mutation.to_dict())
    
    make_parent_directory(out_file)
    
    fh = open(out_file, 'w')
    
    yaml.dump(config, fh)
    
    fh.close()

#=======================================================================================================================
# Post processing code
#=======================================================================================================================
def build_table(args):
    _build_table(
        args.config_file, 
        args.out_file, 
        args.burnin, 
        args.mesh_size, 
        args.table_type, 
        args.thin
    )

def _build_table(config_file, out_file, burnin, mesh_size, table_type, thin):
    if table_type == 'cluster':
        df = post_process.clusters.load_table(
            config_file, 
            burnin=burnin, 
            thin=thin, 
            mesh_size=mesh_size
        )
    
    elif table_type == 'loci':
        df = post_process.loci.load_table(
            config_file, 
            burnin, 
            thin, 
            old_style=False
        )
    
    elif table_type == 'old_style':
        df = post_process.loci.load_table(
            config_file, 
            burnin, 
            thin, 
            old_style=True
        )
        
    df.to_csv(out_file, index=False, sep='\t')

def cluster_plot(args):
    _cluster_plot(
        args.config_file, 
        args.plot_file, 
        args.burnin, 
        args.mesh_size,
        args.min_cluster_size,
        args.plot_type, 
        args.samples,
        args.thin, 
    )

def _cluster_plot(config_file, plot_file, burnin, mesh_size, min_cluster_size, plot_type, samples, thin):
    
    if plot_type == 'density':
        
        plot.clusters.density_plot(
            config_file, 
            plot_file, 
            burnin=burnin, 
            thin=thin, 
            mesh_size=mesh_size,
            min_cluster_size=min_cluster_size,
            samples=samples, 
        )
    
    elif plot_type == 'parallel_coordinates':
        
        plot.clusters.parallel_coordinates_plot(
            config_file, 
            plot_file, 
            burnin=burnin, 
            mesh_size=mesh_size,
            min_cluster_size=min_cluster_size,
            samples=samples, 
            thin=thin
        )
    
    elif plot_type == 'scatter':
        
        plot.clusters.scatter_plot(
            config_file, 
            plot_file, 
            burnin=burnin, 
            mesh_size=mesh_size, 
            min_cluster_size=min_cluster_size, 
            samples=samples, 
            thin=thin
        )
        
def loci_plot(args):
    _loci_plot(
        args.config_file, 
        args.plot_file, 
        args.plot_type, 
        burnin=args.burnin,
        min_cluster_size=args.min_cluster_size,
        samples=args.samples, 
        thin=args.thin)

def _loci_plot(
    config_file, 
    plot_file, 
    plot_type, 
    burnin=0, 
    min_cluster_size=0, 
    samples=None, 
    thin=1):
    
    kwargs = {
        'burnin' : burnin,
        'min_cluster_size' : min_cluster_size,
        'samples' : samples,
        'thin' : thin
    }
    
    if plot_type.startswith('vaf'):
        
        kwargs['value'] = 'variant_allele_frequency'
    
    if plot_type == 'density':
        
        plot.loci.density_plot(
            config_file, 
            plot_file, 
            **kwargs
        )
    
    elif plot_type == 'factor':
        
        plot.loci.factor_plot(
            config_file, 
            plot_file, 
            **kwargs
        )
    
    elif plot_type.endswith('parallel_coordinates'):

        plot.loci.parallel_coordinates_plot(
            config_file, 
            plot_file, 
            **kwargs
        )
    
    elif plot_type.endswith('scatter'):

        plot.loci.scatter_plot(
            config_file, 
            plot_file, 
            **kwargs
        )
    
    elif plot_type == 'similarity_matrix':
        
        plot.loci.similarity_matrix_plot(
            config_file, 
            plot_file, 
            **kwargs
        )

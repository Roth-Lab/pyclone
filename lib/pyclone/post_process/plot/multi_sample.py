'''
Created on 2013-06-06

@author: Andrew Roth
'''
from __future__ import division

import matplotlib.pyplot as pp
import os
import pandas as pd
import seaborn as sb
import yaml

from pyclone.config import load_mutation_from_dict
from pyclone.post_process.cluster import cluster_pyclone_trace
from pyclone.post_process.utils import load_cellular_frequencies_trace

from .utils import setup_axes, setup_plot

def plot_clusters(config_file, plot_file, y_value, burnin=0, thin=1, samples=None, separate_lines=False):
    setup_plot()
    
    if samples is None:
        with open(config_file) as fh:
            samples = yaml.load(fh)['samples'].keys()
    
    df = load_multi_sample_table(config_file, burnin, thin)
    
    if separate_lines:
        plot_df = df
        
    else:
        plot_df = df.groupby(['sample', 'cluster_id']).mean().reset_index()
    
    plot_df['sample_index'] = plot_df['sample'].apply(lambda x: samples.index(x))
    
    plot_df = plot_df.sort_values(by='sample_index')

    xticks = range(len(samples))
    
    if separate_lines:
        fig = pp.figure()
        
        ax = fig.add_subplot(1, 1, 1)
        
        sb.tsplot(
            ax=ax,
            
            data=df,
            condition='cluster_id',
            unit='mutation_id',
            time='sample_index', 
            value=y_value,
            
            err_style='unit_traces',
            marker="o",
            markersize=4
        )
        
    else:
        grid = sb.FacetGrid(
            plot_df,
            hue='cluster_id',
            hue_order=sorted(plot_df['cluster_id'].unique()),
            palette='husl'
        )
        
        ax = grid.ax
        
        fig = grid.fig 
        
        if y_value == 'cellular_prevalence':
            grid.map(
                pp.errorbar, 
                'sample_index', 
                'cellular_prevalence', 
                'cellular_prevalence_std', 
                marker="o",
                markersize=4
            )
        
        elif y_value == 'variant_allele_frequency':
            grid.map(
                pp.plot, 
                'sample_index', 
                'variant_allele_frequency', 
                marker="o",
                markersize=4
            )

    setup_axes(ax)
        
    ax.set_xlim(min(xticks) - 0.1, max(xticks) + 0.1)
    
    ax.set_ylim(0, 1)
    
    ax.set_xlabel('Sample')
    
    if y_value == 'cellular_prevalence':
        ax.set_ylabel('Cellular prevalence')
    
    elif y_value == 'variant_allele_frequency':
        ax.set_ylabel('VAF')
    
    ax.set_xticks(xticks)
    
    ax.set_xticklabels(samples, size=8, rotation=90)
    
    # Shrink current axis by 20%
    box = ax.get_position()
    
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Cluster')

    height = 4
    
    width = 1 * len(samples) + 2
        
    fig.set_figheight(height)
    
    fig.set_figwidth(width)
    
    fig.savefig(plot_file, bbox_inches='tight')

def load_multi_sample_table(config_file, burnin, thin):
    with open(config_file) as fh:
        config = yaml.load(fh)
    
    data = pd.merge(
        _load_variant_allele_frequencies(config),
        _load_cellular_prevalences(config, burnin, thin),
        how='inner',
        on=['mutation_id', 'sample']
    )
    
    trace_dir = os.path.join(config['working_dir'], config['trace_dir'])
    
    labels_file = os.path.join(trace_dir, 'labels.tsv.bz2')
    
    labels = cluster_pyclone_trace(labels_file, burnin, thin)
 
    data = pd.merge(data, labels, on='mutation_id', how='inner')

    return data 

#=======================================================================================================================
# Load allelic prevalences for all samples
#=======================================================================================================================
def _load_variant_allele_frequencies(config):
    data = []
    
    for sample_id in config['samples']:
        file_name = config['samples'][sample_id]['mutations_file']
        
        file_name = os.path.join(config['working_dir'], file_name)
        
        sample_data = _load_sample_variant_allele_frequencies(file_name)
        
        sample_data['sample'] = sample_id
        
        data.append(sample_data)   
    
    data = pd.concat(data, axis=0)
    
    num_samples = len(config['samples'])
    
    # Filter for mutations in all samples
    data = data.groupby('mutation_id').filter(lambda x: len(x) == num_samples)
    
    return data
          
def _load_sample_variant_allele_frequencies(file_name):
    '''
    Load data from PyClone formatted input file.
    '''
    data = {}

    with open(file_name) as fh:
        config = yaml.load(fh)
    
    for mutation_dict in config['mutations']:
        mutation = load_mutation_from_dict(mutation_dict)

        data[mutation.id] = mutation.var_counts / (mutation.ref_counts + mutation.var_counts)

    data = pd.Series(data, name='variant_allele_frequency')
    
    data.index.name = 'mutation_id'
    
    data = data.reset_index()

    return data

#=======================================================================================================================
# Load cellular prevalences for all samples
#=======================================================================================================================
def _load_cellular_prevalences(config, burnin, thin):
    data = []
    
    trace_dir = os.path.join(config['working_dir'], config['trace_dir'])
    
    sample_ids = config['samples'].keys()
    
    for sample_id in sample_ids:
        file_name = os.path.join(trace_dir, '{0}.cellular_frequencies.tsv.bz2'.format(sample_id))
    
        sample_data = _load_sample_cellular_prevalences(file_name, burnin, thin)
    
        sample_data['sample'] = sample_id
    
        data.append(sample_data)
    
    data = pd.concat(data, axis=0)

    num_samples = len(sample_ids)
    
    # Filter for mutations in all samples
    data = data.groupby('mutation_id').filter(lambda x: len(x) == num_samples)
    
    return data 

def _load_sample_cellular_prevalences(file_name, burnin, thin):
    data = load_cellular_frequencies_trace(file_name, burnin, thin)
    
    data = pd.concat([data.mean(), data.std()], axis=1)

    data.columns = 'cellular_prevalence', 'cellular_prevalence_std'
    
    data.index.name = 'mutation_id'
    
    data = data.reset_index()
    
    return data

'''
Created on 2013-06-06

@author: Andrew Roth
'''
from __future__ import division

from collections import OrderedDict
from eppl.parallel_coordinates import aggregated_parallel_coordinates_plot, parallel_coordinates_plot

import csv
import os
import yaml

try:
    import matplotlib.pyplot as plot
except:
    raise Exception("The multi sample plotting module requires the matplotlib package. See http://matplotlib.org.")

try:
    import pandas as pd
except:
    raise Exception("The multi sample plotting module requires the pandas package. See http://http://pandas.pydata.org.")

from pyclone.post_process.cluster import cluster_pyclone_trace
from pyclone.post_process.utils import load_trace

def plot_clusters(config_file, plot_file, prevalence, burnin, thin):
    data = load_results_table(config_file, burnin, thin)
    
    plot_columns = OrderedDict()

    config = _load_yaml_config(config_file)
    
    for entry in config['samples']:
        sample_id = entry['sample_id']
        
        if prevalence == 'allelic':
            col = '{0}_var_freq'.format(sample_id)
        
        elif prevalence == 'cellular':
            col = '{0}_cellular_prevalence_mean'.format(sample_id)
        
        else:
            raise Exception('{0} is an unknown prevalence type.'.format(prevalence))
        
        plot_columns[col] = sample_id

    plot_data = data[plot_columns.keys() + ['cluster_id']]
    
    plot_data.rename(columns=plot_columns, inplace=True)
    
    fig = plot.figure()
    
    ax = fig.add_subplot(1, 1, 1)
    
    title = 'Cluster {0} Prevalence by Sample'.format(prevalence.capitalize())
    
    aggregated_parallel_coordinates_plot(plot_data, 'cluster_id', ax=ax, show_class_size=True,
                                         title=title, x_label='Sample ID', y_label='Prevalence')
    
    ax.set_ylim(0, 1.0)
    
    ax.legend_.set_title('Cluster ID')
    
    fig.savefig(plot_file, dpi=600)
    
    return data

def plot_mutations(config_file, plot_file, prevalence, burnin, thin):
    data = load_results_table(config_file, prevalence, burnin, thin)
    
    fig = plot.figure()
    
    ax = fig.add_subplot(1, 1, 1)
    
    title = 'Mutation {0} Prevalence by Sample'.format(prevalence.capitalize())
    
    parallel_coordinates_plot(data, 'cluster_id', ax=ax, title=title, x_label='Sample ID', y_label='Prevalence')
    
    ax.set_ylim(0, 1.0)
    
    ax.legend_.set_title('Cluster ID')
    
    fig.savefig(plot_file, dpi=600)
    
    return data

def _load_yaml_config(file_name):
    fh = open(file_name)
    
    config = yaml.load(fh)
    
    fh.close()
    
    return config

def load_results_table(config_file, burnin, thin):
    data = []
    
    config = _load_yaml_config(config_file)
    
    trace_dir = os.path.join(config['out_dir'], 'trace')
    
    for entry in config['samples']:
        sample_id = entry['sample_id']
        
        # Allelic prevalence data
        mutations_file = entry['mutations_file']

        sample_allelic_prevalence_data = _load_sample_allelic_prevalences(mutations_file)       
    
        # Cellular prevalence data
        file_name = os.path.join(trace_dir, '{0}.cellular_prevalences.tsv.bz2'.format(sample_id))
        
        sample_cellular_prevalence_data = _load_sample_cellular_prevalences(file_name, burnin, thin)
        
        # Join data frames and rename columns
        sample_data = pd.concat([sample_allelic_prevalence_data, sample_cellular_prevalence_data], axis=1, join='inner')
        
        new_names = ['{0}_{1}'.format(sample_id, x) for x in sample_data.columns]
        
        new_name_map = dict(zip(sample_data.columns, new_names))
        
        sample_data.rename(columns=new_name_map, inplace=True)
        
        data.append(sample_data)
    
    # Join up data from all samples
    data = pd.concat(data, axis=1, join='inner')
    
    # Add cluster cluster labels
    trace_dir = os.path.join(config['out_dir'], 'trace')
    
    labels_file = os.path.join(trace_dir, 'labels.tsv.bz2')
    
    labels = cluster_pyclone_trace(labels_file, burnin, thin)
        
    data = pd.concat([data, labels], axis=1, join='inner')
    
    return data 

#=======================================================================================================================
# Load allelic prevalences for all samples
#=======================================================================================================================
def _load_sample_allelic_prevalences(file_name):
    '''
    Load data from PyClone formatted input file.
    '''
    fh = open(file_name)
    
    data = pd.read_csv(fh, index_col=0, sep='\t')
    
    fh.close()
    
    data['var_freq'] = data['var_counts'] / (data['ref_counts'] + data['var_counts'])
    
    return data

def _load_sample_cellular_prevalences(file_name, burnin, thin):
    trace = load_trace(file_name, burnin, thin)
    
    cellular_prevalence_mean = trace.mean(axis=1)
    
    cellular_prevalence_variance = trace.var(axis=1)

    return pd.DataFrame({'cellular_prevalence_mean' : cellular_prevalence_mean,
                         'cellular_prevalence_variance' : cellular_prevalence_variance}) 

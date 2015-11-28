'''
Created on Nov 27, 2015

@author: Andrew Roth
'''
import os
import pandas as pd
import yaml

from pyclone.config import load_mutation_from_dict

from .cluster import cluster_pyclone_trace
from .utils import load_cellular_frequencies_trace

def load_multi_sample_table(config_file, burnin, thin, old_style=False):
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
    
    if old_style:
        data = _reformat_multi_sample_table(data)

    return data

def _reformat_multi_sample_table(df):
    mean_df = df[['mutation_id', 'sample', 'cellular_prevalence']]
    
    mean_df = mean_df.pivot(index='mutation_id', columns='sample', values='cellular_prevalence')
    
    std_df = df[['mutation_id', 'sample', 'cellular_prevalence_std']]
    
    std_df = std_df.pivot(index='mutation_id', columns='sample', values='cellular_prevalence_std')
    
    std_df = std_df.rename(columns=lambda x: '{0}_std'.format(x))

    cluster_df = df[['mutation_id', 'cluster_id']]
    
    cluster_df = cluster_df.groupby('mutation_id')['cluster_id'].apply(lambda x: x.iloc[0])
    
    return pd.concat([mean_df, std_df, cluster_df], axis=1).reset_index()

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
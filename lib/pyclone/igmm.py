'''
Created on 2013-04-28

@author: Andrew Roth
'''
from __future__ import division

from collections import OrderedDict
from pydp.base_measures import GaussianGammaBaseMeasure
from pydp.data import GaussianData
from pydp.densities import GaussianDensity
from pydp.samplers.atom import GaussianGammaGaussianAtomSampler
from pydp.samplers.dp import DirichletProcessSampler
from pydp.samplers.partition import AuxillaryParameterPartitionSampler

import csv
import os
import yaml

from pyclone.multi_sample import MultiSampleBaseMeasure, MultiSampleDensity, MultiSampleAtomSampler
from pyclone.trace import DiskTrace

def run_igmm_analysis(config_file, trace_dir, num_iters, alpha, alpha_priors):
    data, sample_ids = _load_data(config_file)
    
    sample_atom_samplers = OrderedDict()
    
    sample_base_measures = OrderedDict()
    
    sample_cluster_densities = OrderedDict()
    
    base_measure_params = _load_base_measure_params(config_file)
    
    for sample_id in sample_ids:
        sample_base_measures[sample_id] = GaussianGammaBaseMeasure(base_measure_params['mean'], 
                                                                   base_measure_params['size'], 
                                                                   base_measure_params['alpha'],
                                                                   base_measure_params['beta'])
        
        sample_cluster_densities[sample_id] = GaussianDensity()
        
        sample_atom_samplers[sample_id] = GaussianGammaGaussianAtomSampler(sample_base_measures[sample_id], 
                                                                           sample_cluster_densities[sample_id])    
    
    base_measure = MultiSampleBaseMeasure(sample_base_measures)
    
    cluster_density = MultiSampleDensity(sample_cluster_densities)
    
    atom_sampler = MultiSampleAtomSampler(base_measure, cluster_density, sample_atom_samplers)
    
    partition_sampler = AuxillaryParameterPartitionSampler(base_measure, cluster_density)
    
    sampler = DirichletProcessSampler(atom_sampler, partition_sampler, alpha, alpha_priors)
    
    trace = DiskTrace(trace_dir, sample_ids, data.keys(), {'cellular_frequencies' : 'mean'})
    
    trace.open()
    
    sampler.sample(data.values(), trace, num_iters)
    
    trace.close()

def _load_data(file_name):
    '''
    Load data for all samples.
    
    Args:
        file_name : (str) Path to YAML format configuration file.
    '''
    fh = open(file_name)
    
    config = yaml.load(fh)
    
    fh.close()
    
    sample_data = OrderedDict()
    
    for sample_id in config['samples']:
        file_name = config['samples'][sample_id]['mutations_file']
        
        file_name = os.path.join(config['working_dir'], file_name)
        
        sample_data[sample_id] = _load_sample_data(file_name)        
    
    sample_ids = sample_data.keys()
    
    common_mutations = set.intersection(*[set(x.keys()) for x in sample_data.values()])
    
    data = OrderedDict()
    
    for mutation_id in common_mutations:
        data[mutation_id] = OrderedDict()
        
        for sample_id in sample_ids:
            data[mutation_id][sample_id] = sample_data[sample_id][mutation_id]
    
    return data, sample_ids

def _load_sample_data(file_name):
    data = OrderedDict()
    
    reader = csv.DictReader(open(file_name), delimiter='\t')
    
    for row in reader:
        mutation_id = row['mutation_id']
        
        a = int(row['ref_counts'])
        
        b = int(row['var_counts'])
        
        d = a + b
        
        f = b / d
        
        data[mutation_id] = GaussianData(f)
    
    return data

def _load_base_measure_params(file_name):
    fh = open(file_name)
    
    config = yaml.load(fh)
    
    fh.close()
    
    params = config['base_measure_params']
    
    return params
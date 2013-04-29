'''
Created on 2013-04-28

@author: Andrew Roth
'''
from __future__ import division

from collections import OrderedDict
from pydp.base_measures import BetaBaseMeasure, GammaBaseMeasure
from pydp.data import BinomialData, GammaData
from pydp.densities import Density, log_beta_binomial_pdf
from pydp.samplers.atom import BaseMeasureAtomSampler
from pydp.samplers.dp import DirichletProcessSampler
from pydp.samplers.partition import AuxillaryParameterPartitionSampler

import csv
import os
import yaml

from pydp.samplers.global_params import MetropolisHastingsGlobalParameterSampler
from pydp.proposal_functions import GammaProposal

from pyclone.multi_sample import MultiSampleBaseMeasure, MultiSampleDensity, MultiSampleAtomSampler,\
                                 MultiSampleProposalFunction
from pyclone.trace import DiskTrace

def run_ibbmm_analysis(config_file, trace_dir, num_iters, alpha, alpha_priors):
    data, sample_ids = _load_data(config_file)
    
    sample_atom_samplers = OrderedDict()
    
    sample_base_measures = OrderedDict()
    
    sample_cluster_densities = OrderedDict()
    
    sample_proposal_funcs = OrderedDict()
    
    sample_global_prior = OrderedDict()
    
    base_measure_params = _load_base_measure_params(config_file)
    
    for sample_id in sample_ids:
        sample_base_measures[sample_id] = BetaBaseMeasure(base_measure_params['alpha'], base_measure_params['beta'])
        
        sample_cluster_densities[sample_id] = BetaBinomialDensity(GammaData(100))
        
        sample_atom_samplers[sample_id] = BaseMeasureAtomSampler(sample_base_measures[sample_id], 
                                                                 sample_cluster_densities[sample_id])
        
        sample_proposal_funcs[sample_id] = GammaProposal(1e-2)
        
        sample_global_prior[sample_id] = GammaBaseMeasure(1e-2, 1e-2)
    
    base_measure = MultiSampleBaseMeasure(sample_base_measures)
    
    cluster_density = MultiSampleDensity(sample_cluster_densities, shared_params=True)
    
    proposal_func = MultiSampleProposalFunction(sample_proposal_funcs)
    
    global_prior = MultiSampleBaseMeasure(sample_global_prior)
    
    atom_sampler = MultiSampleAtomSampler(base_measure, cluster_density, sample_atom_samplers)
    
    partition_sampler = AuxillaryParameterPartitionSampler(base_measure, cluster_density)
    
    global_params_sampler = MetropolisHastingsGlobalParameterSampler(GammaBaseMeasure(1e-2, 1e-2), cluster_density, GammaProposal(1e-2))
    
    sampler = DirichletProcessSampler(atom_sampler, partition_sampler, alpha, alpha_priors, global_params_sampler)
    
    trace = DiskTrace(trace_dir, sample_ids, data.keys(), {'cellular_frequencies' : 'x'})
    
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
        
        data[mutation_id] = BinomialData(b, d)
    
    return data 

def _load_base_measure_params(file_name):
    fh = open(file_name)
    
    config = yaml.load(fh)
    
    fh.close()
    
    params = config['base_measure_params']
    
    return params

class BetaBinomialDensity(Density):
    def _log_p(self, data, params):
        a = params.x * self.params.x
        b = self.params.x - a
        
        return log_beta_binomial_pdf(data.x, data.n, a, b)
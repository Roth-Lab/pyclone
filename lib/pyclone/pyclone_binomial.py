'''
Created on 2013-04-28

@author: Andrew Roth
'''
from __future__ import division

from collections import OrderedDict, namedtuple
from math import log
from pydp.base_measures import BetaBaseMeasure
from pydp.densities import log_binomial_pdf, Density
from pydp.samplers.atom import BaseMeasureAtomSampler
from pydp.samplers.dp import DirichletProcessSampler
from pydp.samplers.partition import AuxillaryParameterPartitionSampler

import os
import yaml

from pyclone.config import load_mutation_from_dict
from pyclone.multi_sample import MultiSampleBaseMeasure, MultiSampleDensity, MultiSampleAtomSampler
from pyclone.trace import DiskTrace
from pydp.utils import log_sum_exp

def run_pyclone_binomial_analysis(config_file, trace_dir, num_iters, alpha, alpha_priors):
    data, sample_ids, tumour_content = _load_data(config_file)
    
    sample_atom_samplers = OrderedDict()
    
    sample_base_measures = OrderedDict()
    
    sample_cluster_densities = OrderedDict()
    
    base_measure_params = _load_base_measure_params(config_file)
    
    for sample_id in sample_ids:
        sample_base_measures[sample_id] = BetaBaseMeasure(base_measure_params['alpha'], base_measure_params['beta'])
        
        sample_cluster_densities[sample_id] = PyCloneBinomialDensity(PyCloneBinomialParameter(tumour_content[sample_id]))
        
        sample_atom_samplers[sample_id] = BaseMeasureAtomSampler(sample_base_measures[sample_id], 
                                                                 sample_cluster_densities[sample_id])  
    
    base_measure = MultiSampleBaseMeasure(sample_base_measures)
    
    cluster_density = MultiSampleDensity(sample_cluster_densities)
    
    atom_sampler = MultiSampleAtomSampler(base_measure, cluster_density, sample_atom_samplers)
    
    partition_sampler = AuxillaryParameterPartitionSampler(base_measure, cluster_density)
    
    sampler = DirichletProcessSampler(atom_sampler, partition_sampler, alpha, alpha_priors)
    
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
    
    tumour_content = OrderedDict()
    
    for sample_id in config['samples']:
        file_name = config['samples'][sample_id]['mutations_file']
        
        file_name = os.path.join(config['working_dir'], file_name)
        
        error_rate = config['samples'][sample_id]['error_rate']
        
        sample_data[sample_id] = _load_sample_data(file_name, error_rate)
        
        tumour_content[sample_id] = config['samples'][sample_id]['tumour_content']['value']       
    
    sample_ids = sample_data.keys()
    
    common_mutations = set.intersection(*[set(x.keys()) for x in sample_data.values()])
    
    data = OrderedDict()
    
    for mutation_id in common_mutations:
        data[mutation_id] = OrderedDict()
        
        for sample_id in sample_ids:
            data[mutation_id][sample_id] = sample_data[sample_id][mutation_id]
    
    return data, sample_ids, tumour_content

def _load_sample_data(file_name, error_rate):
    '''
    Load data from PyClone formatted input file.
    '''
    data = OrderedDict()
    
    fh = open(file_name)
    
    config = yaml.load(fh)
    
    fh.close()

    for mutation_dict in config['mutations']:
        mutation = load_mutation_from_dict(mutation_dict)

        data[mutation.id] = _get_pyclone_data(mutation, error_rate)
    
    return data

def _get_pyclone_data(mutation, error_rate):
    a = mutation.ref_counts
    b = mutation.var_counts
    
    d = a + b 
    
    cn_n = tuple([x.cn_n for x in mutation.states])
    cn_r = tuple([x.cn_r for x in mutation.states])
    cn_v = tuple([x.cn_v for x in mutation.states])
    
    mu_n = tuple([x.get_mu_n(error_rate) for x in mutation.states])
    mu_r = tuple([x.get_mu_r(error_rate) for x in mutation.states])
    mu_v = tuple([x.get_mu_v(error_rate) for x in mutation.states])
    
    prior_weights = tuple([x.prior_weight for x in mutation.states])
    
    log_pi = _get_log_pi(prior_weights)
    
    return PyCloneBinomialData(b, d, cn_n, cn_r, cn_v, mu_n, mu_r, mu_v, log_pi)
    
def _get_log_pi(weights):
    pi = [x / sum(weights) for x in weights]
    
    return tuple([log(x) for x in pi])  

def _load_base_measure_params(file_name):
    fh = open(file_name)
    
    config = yaml.load(fh)
    
    fh.close()
    
    params = config['base_measure_params']
    
    return params

PyCloneBinomialData = namedtuple('PyCloneBinomialData',
                                 ['b', 'd', 'cn_n', 'cn_r', 'cn_v', 'mu_n', 'mu_r', 'mu_v', 'log_pi'])

PyCloneBinomialParameter = namedtuple('PyCloneBinomialParameter', 'tumour_content')

class PyCloneBinomialDensity(Density):
    def _log_p(self, data, params):
        ll = []
        
        for cn_n, cn_r, cn_v, mu_n, mu_r, mu_v, log_pi  in zip(data.cn_n, data.cn_r, data.cn_v, data.mu_n, data.mu_r, data.mu_v, data.log_pi):
            temp = log_pi + self._log_binomial_likelihood(data.b,
                                                          data.d,
                                                          cn_n,
                                                          cn_r,
                                                          cn_v,
                                                          mu_n,
                                                          mu_r,
                                                          mu_v,
                                                          params.x)
            
            ll.append(temp)
        
        return log_sum_exp(ll)
    
    def _log_binomial_likelihood(self, b, d, cn_n, cn_r, cn_v, mu_n, mu_r, mu_v, cellular_prevalence):  
        f = cellular_prevalence
        t = self.params.tumour_content
        
        p_n = (1 - t) * cn_n
        p_r = t * (1 - f) * cn_r
        p_v = t * f * cn_v
        
        norm_const = p_n + p_r + p_v
        
        p_n = p_n / norm_const
        p_r = p_r / norm_const
        p_v = p_v / norm_const
        
        mu = p_n * mu_n + p_r * mu_r + p_v * mu_v
        
        return log_binomial_pdf(b, d, mu)

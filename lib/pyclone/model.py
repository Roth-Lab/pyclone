'''
Created on 2013-04-28

@author: Andrew Roth
'''
from __future__ import division

from collections import OrderedDict, namedtuple
from math import log
from pydp.base_measures import BetaBaseMeasure, GammaBaseMeasure
from pydp.data import GammaData
from pydp.densities import Density, log_beta_binomial_pdf
from pydp.proposal_functions import GammaProposal
from pydp.samplers.atom import BaseMeasureAtomSampler
from pydp.samplers.dp import DirichletProcessSampler
from pydp.samplers.global_params import MetropolisHastingsGlobalParameterSampler
from pydp.samplers.partition import AuxillaryParameterPartitionSampler
from pydp.utils import log_sum_exp

import yaml

from pyclone.config import load_mutations_from_csv
from pyclone.multi_sample import MultiSampleBaseMeasure, MultiSampleDensity, MultiSampleAtomSampler
from pyclone.trace import DiskTrace

def run_pyclone_analysis(args, config_file, trace_dir, num_iters, alpha, alpha_priors):
    data, sample_ids = _load_data(config_file, args.error_rate)
    
    sample_atom_samplers = OrderedDict()
    
    sample_base_measures = OrderedDict()
    
    sample_cluster_densities = OrderedDict()
    
    for sample_id in sample_ids:
        sample_base_measures[sample_id] = BetaBaseMeasure(args.base_measure_alpha, args.base_measure_beta)
        
        sample_cluster_densities[sample_id] = SubclonalCopyNumberBetaBinomialDensity(GammaData(args.precision_init))
        
        sample_atom_samplers[sample_id] = BaseMeasureAtomSampler(sample_base_measures[sample_id], 
                                                                 sample_cluster_densities[sample_id])  
    
    base_measure = MultiSampleBaseMeasure(sample_base_measures)
    
    cluster_density = MultiSampleDensity(sample_cluster_densities, shared_params=True)
    
    atom_sampler = MultiSampleAtomSampler(base_measure, cluster_density, sample_atom_samplers)
    
    partition_sampler = AuxillaryParameterPartitionSampler(base_measure, cluster_density)
    
    mh_base_measure = GammaBaseMeasure(args.precision_shape, args.precision_rate)
    
    mh_proposal = GammaProposal(args.precision_proposal_precision)
    
    global_params_sampler = MetropolisHastingsGlobalParameterSampler(mh_base_measure,
                                                                     cluster_density,
                                                                     mh_proposal)
    
    sampler = DirichletProcessSampler(atom_sampler, partition_sampler, alpha, alpha_priors, global_params_sampler)
    
    trace = DiskTrace(trace_dir, sample_ids, data.keys(), {'cellular_prevalences' : 'x'}, precision=True)
    
    trace.open()
    
    sampler.sample(data.values(), trace, num_iters)
    
    trace.close()

def _load_data(file_name, error_rate):
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
    
    for entry in config['samples']:
        sample_id = entry['sample_id']
        
        mutations_file = entry['mutations_file']

        tumour_content = entry['tumour_content']
        
        sample_data[sample_id] = _load_sample_data(mutations_file, error_rate, tumour_content)
    
    sample_ids = sample_data.keys()
    
    common_mutations = set.intersection(*[set(x.keys()) for x in sample_data.values()])
    
    data = OrderedDict()
    
    for mutation_id in common_mutations:
        data[mutation_id] = OrderedDict()
        
        for sample_id in sample_ids:
            data[mutation_id][sample_id] = sample_data[sample_id][mutation_id]
    
    return data, sample_ids

def _load_sample_data(file_name, error_rate, tumour_content):
    '''
    Load data from PyClone formatted input file.
    '''
    data = OrderedDict()
    
    mutations = load_mutations_from_csv(file_name)

    for mutation in mutations:
        data[mutation.id] = _get_pyclone_data(mutation, error_rate, tumour_content)
    
    return data

def _get_pyclone_data(mutation, error_rate, tumour_content):
    a = mutation.ref_counts
    
    b = mutation.var_counts
    
    d = a + b 
    
    return SubclonalCopyNumberData(b, 
                                   d, 
                                   tumour_content, 
                                   mutation.cn_prevalence, 
                                   mutation.normal_cn, 
                                   mutation.major_cn, 
                                   mutation.minor_cn,
                                   mutation.total_cn,
                                   error_rate)

SubclonalCopyNumberData = namedtuple('SubclonalCopyNumberData',
                                     ['b', 'd', 'tumour_content', 'cn_prevalence', 'normal_cn', 'major_cn', 'minor_cn', 'total_cn', 'error_rate'])

class SubclonalCopyNumberBetaBinomialDensity(Density):
    def _log_p(self, data, params):
        prior_weight = -log(3)
        
        ll = [
              prior_weight + self._log_early_event_likelihood(params.x, data, data.major_cn),
              prior_weight + self._log_early_event_likelihood(params.x, data, data.minor_cn),
              prior_weight + self._log_late_event_likelihood(params.x, data)
              ]
        
        return log_sum_exp(ll)
    
    def _log_early_event_likelihood(self, cellular_frequency, data, mutation_cn):
        s = data.cn_prevalence
        
        t = data.tumour_content
        
        f = cellular_frequency
        
        p = [
             (1 - t) * data.normal_cn,
             t * (1 - f) * data.normal_cn,
             t * (1 - s) * f * data.total_cn,
             t * s * f * data.total_cn
             ]
        
        norm_const = sum(p)
        
        p = [p_i / norm_const for p_i in p]
        
        m = [
             data.error_rate,
             data.error_rate,
             1 / data.total_cn,
             mutation_cn / data.total_cn
             ]
        
        mu = sum([x * y for x, y in zip(p, m)])
        
        param_a = mu * self.params.x
        
        param_b = (1 - mu) * self.params.x
        
        return log_beta_binomial_pdf(data.b, data.d, param_a, param_b)   
        
    def _log_late_event_likelihood(self, cellular_frequency, data):
        s = data.cn_prevalence
        
        t = data.tumour_content
        
        f = cellular_frequency
        
        p = [
             (1 - t) * data.normal_cn,
             t * (1 - s) * data.normal_cn,
             t * s * (1 - f) * data.total_cn,
             t * s * f * data.total_cn
             ]
        
        norm_const = sum(p)
        
        p = [p_i / norm_const for p_i in p]
        
        m = [
             data.error_rate,
             data.error_rate,
             data.error_rate,
             1 / data.total_cn
             ]
        
        mu = sum([x * y for x, y in zip(p, m)])
        
        param_a = mu * self.params.x
        
        param_b = (1 - mu) * self.params.x
        
        return log_beta_binomial_pdf(data.b, data.d, param_a, param_b)
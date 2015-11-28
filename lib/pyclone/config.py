'''
Created on 2013-02-12

@author: Andrew Roth
'''
from __future__ import division

from collections import namedtuple, OrderedDict
from math import log

import pyclone.paths as paths

#=======================================================================================================================
# Load data for sampler
#=======================================================================================================================
PyClonelData = namedtuple(
    'PyCloneBinomialData',
    [
        'b', 
        'd', 
        'tumour_conten', 
        'cn_n', 
        'cn_r', 
        'cn_v', 
        'mu_n', 
        'mu_r', 
        'mu_v', 
        'log_pi'
    ]
)

def load_base_measure_params(config_file):
    config = paths.load_config(config_file)    
    
    params = config['base_measure_params']
    
    return params


def load_precision_params(config_file):
    config = paths.load_config(config_file)
        
    return config['beta_binomial_precision_params']    

def load_data(config_file):
    '''
    Load data for all samples.
    
    Args:
        config_file : (str) Path to YAML format configuration file.
    '''
    sample_data = OrderedDict()
    
    error_rate = paths.get_error_rates(config_file)
    
    tumour_content = paths.get_tumour_contents(config_file)
    
    for sample_id, file_name in paths.get_mutations_files(config_file).items():
        sample_data[sample_id] = _load_sample_data(file_name, error_rate[sample_id], tumour_content[sample_id])
    
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
    
    config = paths.load_config(file_name)
  
    for mutation_dict in config['mutations']:
        mutation = load_mutation_from_dict(mutation_dict)

        data[mutation.id] = _get_pyclone_data(mutation, error_rate, tumour_content)
    
    return data

def _get_pyclone_data(mutation, error_rate, tumour_content):
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
    
    return PyClonelData(b, d, tumour_content, cn_n, cn_r, cn_v, mu_n, mu_r, mu_v, log_pi)

def _get_log_pi(weights):
    pi = [x / sum(weights) for x in weights]
    
    return tuple([log(x) for x in pi])

#=======================================================================================================================
# Parse mutation dict
#=======================================================================================================================
def get_mutation(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn, ref_prior, var_prior):
    states = _get_states(normal_cn, minor_cn, major_cn, ref_prior, var_prior)
    
    mutation = Mutation(mutation_id, ref_counts, var_counts)
    
    for (g_n, g_r, g_v) in states:
        state = State(g_n, g_r, g_v, 1)
        
        mutation.add_state(state)
    
    return mutation

def _get_states(normal_cn, minor_cn, major_cn, ref_prior, var_prior):
    if var_prior == 'parental_copy_number':
        states = _get_parental_copy_number_states(normal_cn, minor_cn, major_cn)
     
    elif var_prior == 'total_copy_number':
        total_cn = minor_cn + major_cn
        
        states = _get_total_copy_number_states(normal_cn, total_cn, ref_prior)
     
    elif var_prior == 'no_zygosity':
        total_cn = minor_cn + major_cn
        
        states = _get_no_zygosity_states(normal_cn, total_cn, ref_prior)
     
    elif var_prior == 'AB':
        states = _get_AB_states()
     
    elif var_prior == 'BB':
        states = _get_BB_states()
     
    else:
        raise Exception('{0} is not a recognised method for setting variant priors.'.format(var_prior))
     
    return states

def _get_parental_copy_number_states(normal_cn, minor_cn, major_cn):
    states = set()
    
    g_n = 'A' * normal_cn
    
    total_cn = minor_cn + major_cn
    
    if normal_cn == total_cn:
        if minor_cn == 0:            
            ref_genotypes = [g_n, g_n]
            
            var_genotypes = [
                             'A' * minor_cn + 'B' * major_cn,
                             'A' * (total_cn - 1) + 'B'                             
                             ]
        else:
            ref_genotypes = [g_n, g_n, g_n]
            
            var_genotypes = [
                             'A' * minor_cn + 'B' * major_cn,
                             'A' * major_cn + 'B' * minor_cn,
                             'A' * (total_cn - 1) + 'B'                             
                             ]
    
    else:
        if minor_cn == 0:            
            ref_genotypes = [g_n, g_n, 'A' * total_cn]
            
            var_genotypes = [
                             'A' * minor_cn + 'B' * major_cn,
                             'A' * (total_cn - 1) + 'B',
                             'A' * (total_cn - 1) + 'B'                          
                             ]
        else:
            ref_genotypes = [g_n, g_n, g_n, 'A' * total_cn]
            
            var_genotypes = [
                             'A' * minor_cn + 'B' * major_cn,
                             'A' * major_cn + 'B' * minor_cn,
                             'A' * (total_cn - 1) + 'B',
                             'A' * (total_cn - 1) + 'B'                      
                             ]
    
    for g_r, g_v in zip(ref_genotypes, var_genotypes):
        state = (g_n, g_r, g_v)
        
        states.add(state)
    
    return sorted(states)

def _get_total_copy_number_states(normal_cn, total_cn, ref_prior):
    states = set()
    
    g_n = 'A' * normal_cn
    
    var_genotypes = []
    
    for num_var_alleles in range(1, total_cn + 1):
        g_v = 'A' * (total_cn - num_var_alleles) + 'B' * num_var_alleles
        
        var_genotypes.append(g_v)
    
    if ref_prior == 'normal':
        ref_genotypes = [g_n, ]
    
    elif ref_prior == 'variant':
        ref_genotypes = ['A' * total_cn, ]
    
    elif ref_prior == 'normal_variant':
        ref_genotypes = set([g_n, 'A' * total_cn])
    
    else:
        raise Exception('{0} is not a recognised method for setting reference population priors.'.format(ref_prior))
    
    for g_r in ref_genotypes:
        for g_v in var_genotypes:
            state = (g_n, g_r, g_v)
            
            states.add(state)
    
    return sorted(states)

def _get_no_zygosity_states(normal_cn, total_cn, ref_prior):
    states = set()
    
    g_n = 'A' * normal_cn
    
    if normal_cn == total_cn:
        ref_genotypes = [g_n, ]
        
        var_genotypes = ['A' * (total_cn - 1) + 'B', ]
    
    else:
        if ref_prior == 'normal':
            ref_genotypes = [g_n, ]
            
            var_genotypes = ['A' * (total_cn - 1) + 'B']
        
        elif ref_prior == 'variant':
            ref_genotypes = ['A' * total_cn, ]
            
            var_genotypes = ['A' * (total_cn - 1) + 'B']
        
        elif ref_prior == 'normal_variant':
            ref_genotypes = [g_n, 'A' * total_cn]
            
            var_genotypes = ['A' * (total_cn - 1) + 'B', 'A' * (total_cn - 1) + 'B']
    
        else:
            raise Exception('{0} is not a recognised method for setting reference population priors.'.format(ref_prior))
            
    for g_r, g_v in zip(ref_genotypes, var_genotypes):
        state = (g_n, g_r, g_v)
        
        states.add(state)
    
    return sorted(states)

def _get_AB_states():
    return [('AA', 'AA', 'AB'), ]

def _get_BB_states():
    return [('AA', 'AA', 'BB'), ]

#=======================================================================================================================
# Helper classes
#=======================================================================================================================
class Mutation(object):
    def __init__(self, mutation_id, ref_counts, var_counts):
        self.id = mutation_id
        
        self.ref_counts = ref_counts
        
        self.var_counts = var_counts
        
        self.states = []
        
    @property
    def cn_n(self):
        return [x.cn_n for x in self.states]
        
    @property
    def cn_r(self):
        return [x.cn_r for x in self.states]
    
    @property
    def cn_v(self):
        return [x.cn_v for x in self.states]
    
    @property
    def prior_weights(self):
        return [x.prior_weight for x in self.states]
    
    def add_state(self, state):
        self.states.append(state)
        
    def get_mu_n(self, error_rate):
        return [x.get_mu_n(error_rate) for x in self.states]
    
    def get_mu_r(self, error_rate):
        return [x.get_mu_r(error_rate) for x in self.states]
    
    def get_mu_v(self, error_rate):
        return [x.get_mu_v(error_rate) for x in self.states]        
    
    def to_dict(self):
        return {
                'id' : self.id,
                'ref_counts' : self.ref_counts,
                'var_counts' : self.var_counts,
                'states' : [x.to_dict() for x in self.states]
                }

class State(object):
    def __init__(self, g_n, g_r, g_v, prior_weight):
        self.g_n = g_n
        
        self.g_r = g_r
        
        self.g_v = g_v
        
        self.prior_weight = prior_weight
    
    @property
    def cn_n(self):
        return len(self.g_n)
    
    @property
    def cn_r(self):
        return len(self.g_r)
    
    @property
    def cn_v(self):
        return len(self.g_v)
    
    def get_mu_n(self, error_rate):
        return self._get_variant_allele_probability(self.g_n, error_rate)
    
    def get_mu_r(self, error_rate):
        return self._get_variant_allele_probability(self.g_r, error_rate)
    
    def get_mu_v(self, error_rate):
        return self._get_variant_allele_probability(self.g_v, error_rate)    
            
    def to_dict(self):
        return {'g_n' : self.g_n, 'g_r' : self.g_r, 'g_v' : self.g_v, 'prior_weight' : self.prior_weight}
    
    def _get_copy_number(self, genotype):
        if genotype is None:
            return 0
        else:
            return len(genotype)
        
    def _get_variant_allele_probability(self, genotype, error_rate):
        if genotype is None:
            return error_rate
        
        num_ref_alleles = genotype.count("A")
        num_var_alleles = genotype.count("B")
        
        cn = len(genotype)
        
        if cn != num_ref_alleles + num_var_alleles:
            raise Exception("{0} is not a valid genotype. Only A or B are allowed as alleles.")
        
        if num_ref_alleles == 0:
            return 1 - error_rate
        elif num_var_alleles == 0:
            return error_rate
        else:
            return num_var_alleles / cn        

#=======================================================================================================================
# Factory functions
#=======================================================================================================================
def load_mutation_from_dict(d):
    mutation_id = d['id']
        
    ref_counts = int(d['ref_counts'])
    var_counts = int(d['var_counts'])
    
    mutation = Mutation(mutation_id, ref_counts, var_counts)
    
    for state_dict in d['states']:
        state = load_state_from_dict(state_dict)
        
        mutation.add_state(state)
    
    return mutation
        
def load_state_from_dict(d):
    g_n = d['g_n']
    g_r = d['g_r']
    g_v = d['g_v']
    
    prior_weight = float(d['prior_weight'])
    
    return State(g_n, g_r, g_v, prior_weight)

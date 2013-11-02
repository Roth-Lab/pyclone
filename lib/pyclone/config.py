'''
Created on 2013-02-12

@author: Andrew Roth
'''
from __future__ import division

def get_mutation(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn, ref_prior, var_prior):
    states = _get_states(normal_cn, minor_cn, major_cn, ref_prior, var_prior)
    
    mutation = Mutation(mutation_id, ref_counts, var_counts)
    
    for (g_n, g_r, g_v) in states:
        state = State(g_n, g_r, g_v, 1)
        
        mutation.add_state(state)
    
    return mutation

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
    def __init__(self, g_n, g_r, g_v, cn_prevalence, prior_weight):
        self.g_n = g_n
        
        self.g_r = g_r
        
        self.g_v = g_v
        
        self.cn_prevalence = cn_prevalence
        
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
        return {
                'g_n' : self.g_n, 
                'g_r' : self.g_r, 
                'g_v' : self.g_v,
                'cn_prevalence' : self.cn_prevalence, 
                'prior_weight' : self.prior_weight}
    
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
    
    if 'cn_prevalence' in  d:
        s = float(d['cn_prevalence'])
    
    else:
        s = 1.0        
    
    return State(g_n, g_r, g_v, s, prior_weight)

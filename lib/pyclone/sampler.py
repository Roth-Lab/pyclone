'''
Created on 2011-12-15

@author: Andrew Roth
'''
from __future__ import division

from collections import OrderedDict
from math import log

from pydp.utils import log_sum_exp
from pydp.densities import log_binomial_pdf, Density
from pydp.base_measures import BetaBaseMeasure
from pydp.samplers.atom import BaseMeasureAtomSampler
from pydp.samplers.dp import DirichletProcessSampler
from pydp.samplers.partition import AuxillaryParameterPartitionSampler

import random

class PyCloneSampler(object):
    def __init__(self, alpha, alpha_shape, alpha_rate):
        base_measure = BetaBaseMeasure(1, 1)
        
        cluster_density = PyCloneDensity()
        
        atom_sampler = BaseMeasureAtomSampler(base_measure, cluster_density)
        
        partition_sampler = AuxillaryParameterPartitionSampler(base_measure, cluster_density)
        
        self.sampler = DirichletProcessSampler(atom_sampler, 
                                               partition_sampler, 
                                               alpha=alpha, 
                                               alpha_shape=alpha_shape, 
                                               alpha_rate=alpha_rate)
    
    def sample(self, data, trace, num_iters, seed, print_freq=100):
        random.seed(seed)
        
        self.sampler.sample(data, trace, num_iters, print_freq)

    

class PyCloneData(object):
    def __init__(self, a, b, states, tumour_content, error_rate):
        self.a = a
        self.b = b
        
        self.d = a + b
        
        self.tumour_content = tumour_content
        
        self.cn_n = [x.cn_n for x in states]
        self.cn_r = [x.cn_r for x in states]
        self.cn_v = [x.cn_v for x in states]
        
        self.mu_n = [x.get_mu_n(error_rate) for x in states]
        self.mu_r = [x.get_mu_r(error_rate) for x in states]
        self.mu_v = [x.get_mu_v(error_rate) for x in states]
        
        prior_weights = [x.prior_weight for x in states]
        
        self.log_pi = self._get_log_pi(prior_weights)
        
    def log_p(self, params):
        if params not in self.cache:
            self.cache[params] = self._compute_log_p(params)
            
            if len(self.cache) > self.max_cache_size:
                self.cache.popitem(last=False)
        
        return self.cache[params]
    
    def _get_log_pi(self, weights):
        pi = [x / sum(weights) for x in weights]
        
        return [log(x) for x in pi]  

class PyCloneDensity(Density):
    def __init__(self):
        self.cache = OrderedDict()
        
        self.max_cache_size = 10000
        
    def log_p(self, data, params):
        key = (data, params)
        
        if key not in self.cache:
            self.cache[key] = self._log_p(data, params)
            
            if len(self.cache) > self.max_cache_size:
                self.cache.popitem(last=False)
        
        return self.cache[key]        
    
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
                                                          params.x,
                                                          data.tumour_content)
            
            ll.append(temp)
        
        return log_sum_exp(ll)
    
    def _log_binomial_likelihood(self, b, d, cn_n, cn_r, cn_v, mu_n, mu_r, mu_v, cellular_frequency, tumour_content):
            f = cellular_frequency
            t = tumour_content
            
            p_n = (1 - t) * cn_n
            p_r = t * (1 - f) * cn_r
            p_v = t * f * cn_v
            
            norm_const = p_n + p_r + p_v
            
            p_n = p_n / norm_const
            p_r = p_r / norm_const
            p_v = p_v / norm_const
            
            mu = p_n * mu_n + p_r * mu_r + p_v * mu_v
            
            return log_binomial_pdf(b, d, mu)    

'''
Created on 2011-12-15

@author: Andrew Roth
'''
from __future__ import division

from collections import OrderedDict
from math import log

from pydp.base_measures import BetaBaseMeasure
from pydp.data import BetaData
from pydp.densities import log_binomial_pdf, Density
from pydp.rvs import inverse_sample_rvs
from pydp.samplers.atom import AtomSampler, BaseMeasureAtomSampler
from pydp.samplers.dp import DirichletProcessSampler
from pydp.samplers.partition import AuxillaryParameterPartitionSampler
from pydp.utils import log_sum_exp

import random

class PyCloneSampler(object):
    def __init__(self, alpha, alpha_shape, alpha_rate, cellular_frequency_mesh=None, tumour_content_updater=None):
        base_measure = BetaBaseMeasure(1, 1)
        
        self.cluster_density = PyCloneDensity()
        
        if cellular_frequency_mesh is None:
            atom_sampler = BaseMeasureAtomSampler(base_measure, self.cluster_density)
        else:
            atom_sampler = PyCloneAtomSampler(base_measure, self.cluster_density, cellular_frequency_mesh)
        
        self.tumour_content_updater = tumour_content_updater
        
        partition_sampler = AuxillaryParameterPartitionSampler(base_measure, self.cluster_density)
        
        self.sampler = DirichletProcessSampler(atom_sampler,
                                               partition_sampler,
                                               alpha=alpha,
                                               alpha_shape=alpha_shape,
                                               alpha_rate=alpha_rate)
    
    def sample(self, data, trace, num_iters, seed, print_freq=100):
        random.seed(seed)
        
        self.sampler.initialise_partition(data)
        
        for i in range(num_iters):
            if i % print_freq == 0:
                print self.sampler.num_iters, self.sampler.partition.number_of_cells, self.sampler.alpha
                
                if self.tumour_content_updater is not None:
                    print self.tumour_content_updater.tumour_content
            
            self.sampler.interactive_sample(data)
            
            state = self.sampler.state
                                    
            if self.tumour_content_updater is not None:
                self.tumour_content_updater.update(data, self.sampler.partition, self.cluster_density)
                
                #state['tumour_content'] = self.tumour_content_updater.tumour_content
            
            trace.update(state)
            
            self.sampler.num_iters += 1
        
class PycloneTumourContentUpdater(object):
    def __init__(self, prior_mean, prior_precision, mesh_size):
        a = prior_mean * prior_precision
        
        b = prior_precision - a
        
        self.base_measure = BetaBaseMeasure(1, 1)
        
        self.tumour_content = prior_mean
        
        self.mesh_size = mesh_size

    def update(self, data, partition, cluster_density):
        def sample_log_p(new_value):
            log_p = self.base_measure.log_p(BetaData(new_value))
            
            for cell in partition.cells:
                for item in cell.items:
                    data_point = data[item]
                    
                    old_value = data_point.tumour_content
                    
                    data_point.tumour_content = new_value
                    
                    log_p += cluster_density.log_p(data[item], cell.value)
                    
                    data_point.tumour_content = old_value
            
            return log_p
        
        tumour_content, _ = inverse_sample_rvs(sample_log_p, 0, 1, self.mesh_size)
        
        for data_point in data:
            data_point.tumour_content = tumour_content
        
        self.tumour_content = tumour_content
        
class PyCloneAtomSampler(AtomSampler):
    '''
    Update the partition values using a Gibbs step. 
    
    Requires a Beta base measure and PyClone data.
    '''
    def __init__(self, base_measure, cluster_density, mesh_size):
        AtomSampler.__init__(self, base_measure, cluster_density)
        
        self.mesh_size = mesh_size
    
    def sample(self, data, partition):
        for cell in partition.cells:
            
            def partition_log_p(x):
                params = BetaData(x)
                
                log_p = self.base_measure.log_p(params)
                
                for item in cell.items:
                    log_p += self.cluster_density.log_p(data[item], params)
                
                return log_p
            
            phi = inverse_sample_rvs(partition_log_p, 0, 1, mesh_size=self.mesh_size)
            
            return BetaData(phi[0])

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
    
    def _get_log_pi(self, weights):
        pi = [x / sum(weights) for x in weights]
        
        return [log(x) for x in pi]  

class PyCloneDensity(Density):
    def __init__(self):
        self.cache = OrderedDict()
        
        self.max_cache_size = 10000
        
    def log_p(self, data, params):     
        return self._log_p(data, params)
    
#         key = (data.a, data.b, params)
#         
#         if key not in self.cache:
#             self.cache[key] = self._log_p(data, params)
#             
#             if len(self.cache) > self.max_cache_size:
#                 self.cache.popitem(last=False)
#         
#         return self.cache[key]        
    
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

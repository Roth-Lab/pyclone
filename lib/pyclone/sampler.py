'''
Created on 2011-12-15

@author: Andrew Roth
'''
from __future__ import division

from collections import namedtuple
from math import log

from pydp.base_measures import BaseMeasure

from pydp.partition import Partition

from pydp.samplers.atom import BaseMeasureAtomSampler
from pydp.samplers.concentration import GammaPriorConcentrationSampler
from pydp.samplers.partition import AuxillaryParameterPartitionSampler

from pydp.rvs import uniform_rvs
from pydp.utils import log_sum_exp

class DirichletProcessSampler(object):
    def __init__(self, cluster_density, tumour_content, alpha=None, alpha_shape=None, alpha_rate=None):
        self.cluster_density = cluster_density
        
        self.base_measure = PyCloneBaseMeasure(tumour_content)

        self.partition_sampler = AuxillaryParameterPartitionSampler(self.base_measure, cluster_density)
        
        self.atom_sampler = BaseMeasureAtomSampler(self.base_measure, cluster_density)           
        
        if alpha is None:
            self.alpha = 1
            
            self.update_alpha = True
            
            self.concentration_sampler = GammaPriorConcentrationSampler(alpha_shape, alpha_rate)
        else:
            self.alpha = alpha
            
            self.update_alpha = False 
        
        self.num_iters = 0
    
    @property
    def state(self):
        return {
                'alpha' : self.alpha,
                'cellular_frequencies' : [param.phi for param in self.partition.item_values],
                'labels' : self.partition.labels,
                'phi' : [param.phi for param in self.partition.cell_values]
                }
    
    def initialise_partition(self, data):
        self.partition = Partition()
        
        for item, _ in enumerate(data):
            self.partition.add_cell(self.base_measure.random())
            
            self.partition.add_item(item, item)        
    
    def sample(self, data, results_db, num_iters, print_freq=1):
        self.initialise_partition(data)
        
        print "Tumour Content :", self.base_measure.tumour_content
        
        for i in range(num_iters):
            if i % print_freq == 0:
                print self.num_iters, self.partition.number_of_cells, self.alpha, self.base_measure.tumour_content
            
            self.interactive_sample(data)
            
            results_db.update_trace(self.state)
            
            self.num_iters += 1
    
    def interactive_sample(self, data):
        if self.update_alpha:
            self.alpha = self.concentration_sampler.sample(self.alpha,
                                                           self.partition.number_of_cells,
                                                           self.partition.number_of_items)
        
        self.partition_sampler.sample(data, self.partition, self.alpha)
        
        self.atom_sampler.sample(data, self.partition)
        
        self._update_tumour_content(data)
    
    def _init_partition(self, base_measure):
        self.partition = Partition()
        
        for item, _ in enumerate(self.data):
            self.partition.add_cell(base_measure.random())
            
            self.partition.add_item(item, item)
            
    def _update_tumour_content(self, data):
        log_f = lambda x: sum([self.cluster_density.log_p(data_point, PyCloneParameter(params.phi, x)) for data_point, params in zip(data, self.partition.item_values)])
        
        self.base_measure.tumour_content, _ = inverse_sample(log_f, 0, 1, 1000)

class PyCloneBaseMeasure(BaseMeasure):
    def __init__(self, tumour_content):
        self.tumour_content = tumour_content
        
    def random(self):
        phi = uniform_rvs(0, 1)
        
        return PyCloneParameter(phi, self.tumour_content)
    
PyCloneParameter = namedtuple('PyCloneParameter', ['phi', 's'])
    
def inverse_sample(log_f, a, b, mesh_size=100):
    u = uniform_rvs(0, 1)
    
    log_u = log(u)

    step_size = (b - a) / mesh_size

    log_step_size = log(b - a) - log(mesh_size)

    knots = [i * step_size + a for i in range(0, mesh_size + 1)]
    
    log_likelihood = [log_f(x) for x in knots]
    
    log_riemann_sum = []
    
    for y in log_likelihood:
        log_riemann_sum.append(y + log_step_size)
    
    log_norm_const = log_sum_exp(log_riemann_sum)
    
    log_cdf = None
    
    for x, y in zip(knots, log_likelihood):
        log_q = y - log_norm_const
        
        log_partial_pdf_riemann_sum = log_q + log_step_size
        
        if log_cdf is None:
            log_cdf = log_partial_pdf_riemann_sum
        else:
            log_cdf = log_sum_exp([log_cdf, log_partial_pdf_riemann_sum])
     
        if log_u < log_cdf:
            break

    return x, log_q    

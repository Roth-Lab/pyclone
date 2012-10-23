'''
Created on 2011-12-15

@author: Andrew Roth
'''
from __future__ import division

from collections import namedtuple
from math import log

from pydp.base_measures import BetaBaseMeasure

from pydp.densities import Density, log_binomial_pdf
from pydp.partition import Partition
from pydp.rvs import uniform_rvs

from pydp.proposal_functions import ProposalFunction

from pydp.samplers.atom import AtomSampler, BaseMeasureAtomSampler, MetropolisHastingsAtomSampler
from pydp.samplers.concentration import GammaPriorConcentrationSampler
from pydp.samplers.partition import AuxillaryParameterPartitionSampler

from pydp.utils import log_sum_exp, memoized
from pydp.data import BetaParameter, BetaData

class DirichletProcessSampler(object):
    def __init__(self, alpha=None):
        self.base_measure = BetaBaseMeasure(1, 1)
        
        cluster_density = PyCloneDensity()
        
        self.partition_sampler = AuxillaryParameterPartitionSampler(self.base_measure, cluster_density)
        
        self.atom_sampler = PyCloneGibbsAtomSampler(self.base_measure, cluster_density)           
        
        if alpha is None:
            self.alpha = 1
            
            self.update_alpha = True
            
            self.concentration_sampler = GammaPriorConcentrationSampler(1, 1)
        else:
            self.alpha = alpha
            
            self.update_alpha = False 
        
        self.num_iters = 0
    
    @property
    def state(self):
        return {
                'alpha' : self.alpha,
                'cellular_frequencies' : [param.x for param in self.partition.item_values],
                'labels' : self.partition.labels,
                'phi' : [param.x for param in self.partition.cell_values]
                }
    
    def initialise_partition(self, data):
        self.partition = Partition()
        
        for item, _ in enumerate(data):
            self.partition.add_cell(self.base_measure.random())
            
            self.partition.add_item(item, item)        
    
    def sample(self, data, results_db, num_iters, print_freq=100):
        self.initialise_partition(data)
        
        for i in range(num_iters):
            if i % print_freq == 0:
                print self.num_iters, self.partition.number_of_cells, self.alpha
            
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
    
    def _init_partition(self, base_measure):
        self.partition = Partition()
        
        for item, _ in enumerate(self.data):
            self.partition.add_cell(base_measure.random())
            
            self.partition.add_item(item, item)

class StromalSampler(DirichletProcessSampler):
    def __init__(self, alpha=None):
        self.base_measure = BetaBaseMeasure(1, 1)
        
        cluster_density = PyCloneStromalDensity()
        
        self.partition_sampler = AuxillaryParameterPartitionSampler(self.base_measure, cluster_density)
        
        proposal_func = PyCloneProposalFunction()
        
        self.atom_sampler = MetropolisHastingsAtomSampler(self.base_measure, cluster_density, proposal_func)
        
        self.s = 0.01           
        
        if alpha is None:
            self.alpha = 1
            
            self.update_alpha = True
            
            self.concentration_sampler = GammaPriorConcentrationSampler(1, 1)
        else:
            self.alpha = alpha
            
            self.update_alpha = False 
        
        self.num_iters = 0
    
    @property
    def state(self):
        return {
                'alpha' : self.alpha,
                'cellular_frequencies' : [param.x for param in self.partition.item_values],
                'labels' : self.partition.labels,
                'phi' : [param.x for param in self.partition.cell_values],
                's' : self.s
                }
        
    def interactive_sample(self, data):
        if self.update_alpha:
            self.alpha = self.concentration_sampler.sample(self.alpha,
                                                           self.partition.number_of_cells,
                                                           self.partition.number_of_items)
        
        self.partition_sampler.sample(data, self.partition, self.alpha)
        
        self.atom_sampler.sample(data, self.partition)
        
        self.update_normal_proportion(data, self.partition)
        
    def update_normal_proportion(self, data, partition):
        old_ll = 0
        new_ll = 0
        
        new_s = uniform_rvs(0, 1)
        
        for cell in partition.cells:
            old_param = cell.value
            
            new_param = old_param.copy()
            new_param.s = new_s

            for j in cell.items:
                old_ll += self.cluster_density.log_p(data[j], old_param)
                new_ll += self.cluster_density.log_p(data[j], new_param)
            
        log_numerator = new_ll
        log_denominator = old_ll
            
        log_ratio = log_numerator - log_denominator
            
        u = uniform_rvs(0, 1)
            
        if log_ratio >= log(u):
            self.s = new_s
            
            for cell in partition.cells:
                cell.value.s = new_s
                
    def _init_partition(self, base_measure):
        self.partition = Partition()
        
        for item, _ in enumerate(self.data):
            self.partition.add_cell(base_measure.random())
            
            self.partition.add_item(item, item)                   

class PyCloneProposalFunction(ProposalFunction):
    def log_p(self, data, params):
        return 0
    
    def random(self, params):
        phi = uniform_rvs(0, 1)
        
        return PyCloneParameter(phi, params.s)

class PyCloneGibbsAtomSampler(AtomSampler):
    def sample(self, data, partition):
        for cell in partition.cells:
            log_f = lambda x: sum([self.cluster_density.log_p(data[item], x) for item in cell.items])
            
            phi = self.inverse_sample(log_f, 0, 1, 100)
            
            cell.value = phi
            
    def inverse_sample(self, log_f, a, b, mesh_size):
        u = uniform_rvs(0, 1)
        
        log_u = log(u)
        
        step_size = (b - a) / mesh_size
        
        log_step_size = log(b - a) - log(mesh_size)
    
        knots = [i * step_size + a for i in range(0, mesh_size + 1)]
        knots = [BetaData(x) for x in knots]
    
        log_total = [log_f(x) + log_step_size for x in knots]
    
        log_norm_const = log_sum_exp(log_total)
        
        log_cdf = None
        
        for x in knots:
            term = log_f(x) - log_norm_const
            
            if log_cdf is None:
                log_cdf = term
            else:
                log_cdf = log_sum_exp([log_cdf, term + log_step_size])

            if log_u < log_cdf:
                return x
        
        return x            

#=======================================================================================================================
# Data class
#=======================================================================================================================
PyCloneData = namedtuple('PyCloneData', ['a', 'd', 'mu_r', 'mu_v', 'log_pi_r', 'log_pi_v'])
PyCloneParameter = namedtuple('PyCloneParameter', ['phi', 's'])

class PyCloneDensity(Density):
    def log_p(self, data, params):
        ll = []
        
        for mu_r, log_pi_r in zip(data.mu_r, data.log_pi_r):
            for mu_v, log_pi_v in zip(data.mu_v, data.log_pi_v):
                temp = log_pi_r + log_pi_v + self._log_binomial_likelihood(data.a, data.d, params.x, mu_r, mu_v)
                
                ll.append(temp)
        
        return log_sum_exp(ll)
    
    def _log_binomial_likelihood(self, a, d, phi, mu_r, mu_v):
        mu = (1 - phi) * mu_r + phi * mu_v
        
        return log_binomial_pdf(a, d, mu)
    
class PyCloneStromalDensity(Density):
    def log_p(self, data, params):
        ll = []
        
        for mu_r, log_pi_r in zip(data.mu_r, data.log_pi_r):
            for mu_v, log_pi_v in zip(data.mu_v, data.log_pi_v):
                temp = log_pi_r + log_pi_v + self._log_binomial_likelihood(data.a, data.d, params.phi, params.s, mu_r, mu_v)
                
                ll.append(temp)
        
        return log_sum_exp(ll)
    
    def _log_binomial_likelihood(self, a, d, phi, s, mu_r, mu_v):
        mu = (1 - phi * (1 - s)) * mu_r + phi * (1 - s) * mu_v
        
        return log_binomial_pdf(a, d, mu)    

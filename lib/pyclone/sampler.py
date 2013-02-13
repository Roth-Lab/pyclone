'''
Created on 2011-12-15

@author: Andrew Roth
'''
from __future__ import division

from collections import OrderedDict, namedtuple
from math import log, lgamma as log_gamma
from numpy.random import multinomial
from random import betavariate as beta_rvs, gammavariate as gamma_rvs, shuffle, uniform as uniform_rvs

import numpy as np
import random

class DirichletProcessSampler(object):
    def __init__(self, tumour_content, alpha=None, alpha_shape=None, alpha_rate=None):
        self.base_measure = BaseMeasure(tumour_content)

        self.partition_sampler = PartitionSampler(self.base_measure)
        
        self.atom_sampler = AtomSampler(self.base_measure)           
        
        if alpha is None:
            self.alpha = 1
            
            self.update_alpha = True
            
            self.concentration_sampler = ConcentrationSampler(alpha_shape, alpha_rate)
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
    
    def sample(self, data, results_db, num_iters, print_freq=100, seed=None):
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)
        
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
    
    def _init_partition(self, base_measure):
        self.partition = Partition()
        
        for item, _ in enumerate(self.data):
            self.partition.add_cell(base_measure.random())
            
            self.partition.add_item(item, item)
        
#=======================================================================================================================
# DP Samplers
#=======================================================================================================================
class AtomSampler(object):
    '''
    Update the atom values using a Metropolis-Hastings steps with the base measure as a proposal density.
    '''
    def __init__(self, base_measure):
        self.base_measure = base_measure
        
        self.proposal_func = BaseMeasureProposalFunction(base_measure)
        
    def sample(self, data, partition):
        for cell in partition.cells:
            old_param = cell.value
            new_param = self.proposal_func.random(old_param)            

            old_ll = 0
            new_ll = 0

            for j in cell.items:
                old_ll += data[j].log_p(old_param)
                new_ll += data[j].log_p(new_param)
            
            log_numerator = new_ll + self.proposal_func.log_p(old_param, new_param)
            log_denominator = old_ll + self.proposal_func.log_p(new_param, old_param)
            
            log_ratio = log_numerator - log_denominator
            
            u = uniform_rvs(0, 1)
            
            if log_ratio >= log(u):
                cell.value = new_param
                
class ConcentrationSampler(object):
    '''
    Sampler for updating concentration parameter. Uses Gibbs update assuming a gamma prior on the concentration
    parameter.
    '''
    def __init__(self, a, b):
        '''
        Args :
            a : (float) Location parameter of the gamma prior.
            b : (float) Scale parameter of the gamma prior.
        '''
        self.a = a
        self.b = b
    
    def sample(self, old_value, num_clusters, num_data_points):
        a = self.a
        b = self.b
        
        k = num_clusters
        n = num_data_points
        
        eta = beta_rvs(old_value + 1, n)
    
        x = (a + k - 1) / (n * (b - log(eta)))
        
        pi = x / (1 + x)
    
        label = discrete_rvs([pi, 1 - pi])
        
        scale = 1 / (b - log(eta))
                
        if label == 0:
            new_value = gamma_rvs(a + k, scale)
        else:
            new_value = gamma_rvs(a + k - 1, scale)
        
        return new_value                

class PartitionSampler(object):
    '''
    Sampler for updating parition of the data. Uses algorithm 8 of Neal "Sampling Methods For Dirichlet Process Mixture
    Models".
    '''
    
    def __init__(self, base_measure):
        self.base_measure = base_measure
    
    def sample(self, data, partition, alpha, m=2):
        '''
        Sample a new partition. 
        '''
        items = range(len(data))
        
        shuffle(items)
        
        for item in items:
            data_point = data[item]
            
            old_cell_index = partition.labels[item]
            
            partition.remove_item(item, old_cell_index)
            
            if partition.counts[old_cell_index] == 0:
                num_new_tables = m - 1
            else:
                num_new_tables = m
            
            for _ in range(num_new_tables):
                partition.add_cell(self.base_measure.random())
            
            log_p = np.zeros(len(partition))
            
            for i, cell in enumerate(partition.cells):
                cluster_log_p = data_point.log_p(cell.value)
                
                counts = cell.size
                
                if counts == 0:
                    counts = alpha / m
                
                log_p[i] = log(counts) + cluster_log_p
    
            log_p = log_space_normalise(log_p)
            
            p = np.exp(log_p)
            
            if p.sum() < 0.9999:
                print p.sum()
            
            p = p / p.sum()
            
            new_cell_index = discrete_rvs(p)
            
            partition.add_item(item, new_cell_index)
            
            partition.remove_empty_cells()

#=======================================================================================================================
# Auxillary DP classes
#=======================================================================================================================
Parameter = namedtuple('Parameter', ['phi', 's'])

class BaseMeasure(object):
    '''
    Class to sampler from the PyClone base measure.
    '''
    def __init__(self, tumour_content):
        self.tumour_content = tumour_content
        
    def random(self):
        phi = uniform_rvs(0, 1)
        
        return Parameter(phi, self.tumour_content)
    
class BaseMeasureProposalFunction(object):
    '''
    Proposal function for atom sampler which draw values from a base measure.
    '''
    def __init__(self, base_measure):
        self.base_measure = base_measure
        
    def log_p(self, data, params):
        return 0
    
    def random(self, params):
        return self.base_measure.random()
    
class DataPoint(object):
    def __init__(self, a, b, cn_n, cn_r, cn_v, mu_n, mu_v, mu_r, prior_weights):
        self.a = a
        self.b = b
        
        self.d = a + b
        
        self.cn_n = np.array(cn_n)
        self.cn_r = np.array(cn_r)
        self.cn_v = np.array(cn_v)
        
        self.mu_n = np.array(mu_n)
        self.mu_r = np.array(mu_r)
        self.mu_v = np.array(mu_v)
        
        self.log_pi = self._get_log_pi(prior_weights)
        
        self.log_norm_const = log_binomial_coefficient(self.d, self.b)
        
        self.cache = OrderedDict()
        
        self.max_cache_size = 10000
    
    def log_p(self, params):
        if params not in self.cache:
            self.cache[params] = self._compute_log_(params)
            
            if len(self.cache) > self.max_cache_size:
                self.cache.popitem(last=False)
        
        return self.cache[params]
        
    def _compute_log_(self, params):
        phi = params.phi
        s = params.s
        
        b = self.b
        d = self.d

        cn_n = self.cn_n        
        cn_r = self.cn_r
        cn_v = self.cn_v
        
        mu_n = self.mu_n
        mu_r = self.mu_r
        mu_v = self.mu_v
        
        p_n = (1 - s) * cn_n
        p_r = s * (1 - phi) * cn_r
        p_v = s * phi * cn_v
        
        norm_const = p_n + p_r + p_v
        
        p_n = p_n / norm_const
        p_r = p_r / norm_const
        p_v = p_v / norm_const
        
        mu = p_n * mu_n + p_r * mu_r + p_v * mu_v
        
        ll = log_binomial_likelihood(b, d, mu) - self.log_norm_const
        
        return log_sum_exp(ll)
    
    def _get_log_pi(self, weights):
        weights = np.array(weights)
        
        pi = weights / weights.sum()
        
        return np.log(pi)  

#=======================================================================================================================
# Partition Data Structure
#=======================================================================================================================
class Partition(object):
    def __init__(self):
        self.cells = []
    
    def __len__(self):
        return len(self.cells)
    
    @property
    def cell_values(self):
        return [cell.value for cell in self.cells]
    
    @property
    def counts(self):
        return [cell.size for cell in self.cells]
    
    @property
    def item_values(self):
        cell_values = self.cell_values
        labels = self.labels
        
        return [cell_values[i] for i in labels]
    
    @property
    def labels(self):
        labels = [None] * self.number_of_items
        
        for cell_index, cell in enumerate(self.cells):
            for item in cell.items:
                labels[item] = cell_index
        
        return labels
    
    @property
    def number_of_cells(self):
        return len(self.cells)
    
    @property
    def number_of_items(self):
        return sum(self.counts)
        
    def add_cell(self, value):
        self.cells.append(PartitionCell(value))
    
    def add_item(self, item, cell_index):
        self.cells[cell_index].add_item(item)
        
    def get_cell_by_value(self, value):
        for cell in self.cells:
            if cell.value == value:
                return cell
    
    def remove_item(self, item, cell_index):
        self.cells[cell_index].remove_item(item)
        
    def remove_empty_cells(self):
        for cell in self.cells[:]:
            if cell.empty:
                self.cells.remove(cell)

    def copy(self):
        partition = Partition()
    
        for cell_index, cell in enumerate(self.cells):
            partition.add_cell(cell.value)
        
            for item in cell.items:
                partition.add_item(item, cell_index)
        
        return partition

class PartitionCell(object):
    def __init__(self, value):
        self.value = value
        
        self._items = []
    
    @property
    def empty(self):
        if self.size == 0:
            return True
        else:
            return False
    
    @property
    def items(self):
        return self._items[:]
    
    @property
    def size(self):
        return len(self._items)
    
    def add_item(self, item):
        self._items.append(item)
    
    def remove_item(self, item):
        self._items.remove(item)
    
    def __contains__(self, x):
        if x in self._items:
            return True
        else:
            return False

#=======================================================================================================================
# Utility functions
#=======================================================================================================================
def discrete_rvs(p):
    '''
    Return a discrete random variable.
    
    Args:
        p : An iterable containing the probabilites for each class.
        
    Returns:
        x : (int) The index of the class drawn.
    '''
    x = multinomial(1, p)
    
    return x.argmax()

def log_binomial_likelihood(x, n, p):
    return x * np.log(p) + (n - x) * np.log(1 - p)

def log_binomial_coefficient(n, x):
    return log_factorial(n) - log_factorial(n - x) - log_factorial(x)

def log_factorial(n):
    return log_gamma(n + 1)

def log_space_normalise(x):
    log_norm_const = log_sum_exp(x)
    
    return x - log_norm_const

def log_sum_exp(x):
    return np.logaddexp.reduce(x)

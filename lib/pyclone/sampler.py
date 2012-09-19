'''
Created on 2011-12-15

@author: Andrew Roth
'''
from __future__ import division

from math import exp, log, lgamma as log_gamma

from random import betavariate as beta_rvs, gammavariate as gamma_rvs, uniform as uniform_rvs, shuffle

from pyclone.utils import discrete_rvs, log_binomial_likelihood, log_space_normalise, log_sum_exp, memoized

class DirichletProcessSampler(object):
    def __init__(self, data, alpha=None):
        self.data = data
        
        if alpha is None:
            self.alpha = 1
            
            self.update_alpha = True
        else:
            self.alpha = alpha
            
            self.update_alpha = False 
        
        self.num_iters = 0
        
        self._init_partition()
    
    @property
    def state(self):
        return {
                'alpha' : self.alpha,
                'labels' : self.partition.labels,
                'phi' : self.partition.item_values
                }
    
    def sample(self, results_db, num_iters, print_freq=100):
        for i in range(num_iters):
            if i % print_freq == 0:
                print self.num_iters, self.partition.number_of_cells, self.alpha
            
            self.interactive_sample()
            
            results_db.update_trace(self.state)
            
            self.num_iters += 1
    
    def interactive_sample(self):
        self.partition = sample_partition(self.alpha,
                                          self.data,
                                          self.partition,
                                          uniform_proposal)
        
        sample_cell_values(self.data, self.partition, uniform_proposal)
        
        if self.update_alpha:
            self.alpha = sample_alpha(self.alpha,
                                      self.partition.number_of_cells,
                                      self.partition.number_of_items,
                                      1,
                                      1)
    
    def _init_partition(self):
        partition = Partition()
        
        for item, _ in enumerate(self.data):
            partition.add_cell(uniform_rvs(0, 1))
            
            partition.add_item(item, item)
            
        self.partition = partition
#=======================================================================================================================
# Data class
#=======================================================================================================================
class PyCloneData(object):
    def __init__(self, a, d, ref_priors, var_priors):
        self.a = a
        self.d = d
        
        self.mu_r = ref_priors.keys()
        self.mu_v = var_priors.keys()
        
        self.log_pi_r = self._get_log_mix_weights(ref_priors.values())
        self.log_pi_v = self._get_log_mix_weights(var_priors.values())
        
    @memoized
    def log_p(self, phi=None):
        ll = []
        
        for mu_r, log_pi_r in zip(self.mu_r, self.log_pi_r):
            for mu_v, log_pi_v in zip(self.mu_v, self.log_pi_v):
                temp = log_pi_r + log_pi_v + self._log_binomial_likelihood(phi, mu_r, mu_v)
                
                ll.append(temp)
        
        return log_sum_exp(ll)

    def _get_log_mix_weights(self, delta):
        log_denominator = log_gamma(sum(delta) + 1)
        
        log_mix_weights = []
        
        for i, d_i in enumerate(delta):
            log_numerator = log_gamma(d_i + 1)
            
            for j, d_j in enumerate(delta):
                if i != j:
                    log_numerator += log_gamma(d_j)
            
            log_mix_weights.append(log_numerator - log_denominator)
        
        return log_mix_weights
    
    def _log_binomial_likelihood(self, phi, mu_r, mu_v):
        mu = (1 - phi) * mu_r + phi * mu_v
        
        return log_binomial_likelihood(self.a, self.d, mu)

#=======================================================================================================================
# DP Data Structures
#=======================================================================================================================
class Partition(object):
    def __init__(self):
        self.cells = []
    
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
        labels = {}
        
        for cell_index, cell in enumerate(self.cells):
            for j in cell.items:
                labels[j] = cell_index
        
        return [labels[j] for j in sorted(labels)]
    
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
# Sampling Fucntions
#=======================================================================================================================
def uniform_proposal(phi=None):
    return uniform_rvs(0, 1)

def sample_alpha(old_value, k, n, a=1, b=1):
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

def sample_partition(alpha, data, old_partition, proposal_func, m=2):
    '''
    Sample a new partition according to algorithm 8 of Neal "Sampling Methods For Dirichlet Process Mixture Models"
    '''
    partition = old_partition.copy()
    
    items = range(len(data))
    
    shuffle(items)
    
    for item in items:
        data_point = data[item]
        
        old_cluster_label = partition.labels[item]
        
        partition.remove_item(item, old_cluster_label)
        
        if partition.counts[old_cluster_label] == 0:
            num_new_tables = m - 1
        else:
            num_new_tables = m
        
        for _ in range(num_new_tables):
            partition.add_cell(proposal_func())
        
        log_p = []
        
        for cell in partition.cells:
            cluster_log_p = data_point.log_p(cell.value)
            
            counts = cell.size
            
            if counts == 0:
                counts = alpha / m
            
            log_p.append(log(counts) + cluster_log_p)

        log_p = log_space_normalise(log_p)
        
        p = [exp(x) for x in log_p]
        
        new_cluster_label = discrete_rvs(p)
        
        partition.add_item(item, new_cluster_label)
        
        partition.remove_empty_cells()
    
    return partition

#def sample_partition(alpha, data, old_partition, proposal_func):
#    '''
#    Sample a new partition according to algorithm 7 of Neal "Sampling Methods For Dirichlet Process Mixture Models"
#    '''
#    partition = old_partition.copy()
#    
#    n = partition.number_of_items
#    
#    for item, data_point in enumerate(data):
#        old_cluster_label = partition.labels[item]
#        old_value = partition.item_values[item]
#        
#        partition.remove_item(item, old_cluster_label)
#        
#        if partition.counts[old_cluster_label] == 0:
#            p = [x / (n - 1) for x in partition.counts]
#            
#            new_cluster_label = discrete_rvs(p)
#            
#            new_value = partition.cell_values[new_cluster_label]
#            
#            old_ll = data_point.log_p(old_value)
#            new_ll = data_point.log_p(new_value)
#            
#            log_ratio = log(n - 1) - log(alpha) + new_ll - old_ll
#            
#            u = uniform_rvs(0, 1)
#            
#            if log_ratio >= log(u):
#                partition.add_item(item, new_cluster_label)
#            else:
#                partition.add_item(item, old_cluster_label)
#        
#        else:
#            new_value = proposal_func()
#            
#            old_ll = data_point.log_p(old_value)
#            new_ll = data_point.log_p(new_value)
#            
#            log_ratio = log(alpha) - log(n - 1) + new_ll - old_ll
#            
#            u = uniform_rvs(0, 1)
#            
#            if log_ratio >= log(u):
#                partition.add_cell(new_value)
#                
#                cell = partition.get_cell_by_value(new_value)
#                
#                cell.add_item(item)
#            else:
#                partition.add_item(item, old_cluster_label)
#    
#    partition.remove_empty_cells()
#    
#    for item, data_point in enumerate(data):
#        old_cluster_label = partition.labels[item]
#        
#        if partition.cells[old_cluster_label].size == 1:
#            continue
#        
#        partition.remove_item(item, old_cluster_label)
#        
#        log_p = []
#        
#        for cell in partition.cells:
#            cluster_log_p = data_point.log_p(cell.value)
#            
#            counts = cell.size
#            
#            log_p.append(log(counts) + cluster_log_p)
#
#        log_p = log_space_normalise(log_p)
#        
#        p = [exp(x) for x in log_p]
#        
#        new_cluster_label = discrete_rvs(p)
#        
#        partition.add_item(item, new_cluster_label)
#    
#    partition.remove_empty_cells()
#    
#    return partition
        
def sample_cell_values(data, partition, proposal_func):
    for cell in partition.cells:
        _sample_cell_value(data, cell, proposal_func)
      
def _sample_cell_value(data, cell, proposal_func):          
    for _ in range(100):
        old_ll = 0
        new_ll = 0
        
        old_param = cell.value
        new_param = proposal_func(old_param)
        
        for j in cell.items:
            old_ll += data[j].log_p(old_param)
            new_ll += data[j].log_p(new_param)
        
        log_ratio = new_ll - old_ll
        
        u = uniform_rvs(0, 1)
        
        if log_ratio >= log(u):
            cell.value = new_param
            break

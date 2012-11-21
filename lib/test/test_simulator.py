'''
Created on 2012-09-15

@author: andrew
'''
from __future__ import division

import random

from random import betavariate as beta_rvs, gammavariate as gamma_rvs, normalvariate as normal_rvs

from pyclone.utils import discrete_rvs, dirichlet_rvs, binomial_rvs

from pyclone.sampler import PyCloneData, Partition, sample_alpha, sample_partition, sample_cell_values

from pyclone.test.crp import sample_from_crp
from pyclone.test.geweke import compare_trace

def normal_proposal(old_phi):
    while True:
        new_phi = normal_rvs(old_phi, 1)
        
        if 0 <= new_phi <= 1:
            return new_phi

def beta_proposal(m):
    if m == 0 or m == 1:
        return random.random()
    
    s = 10
    
    a = s * m
    b = s - a
    
    if a <= 0:
        a = 1e-100
    if b <= 0:
        b = 1e-100

    return beta_rvs(a, b)  

def uniform_proposal(phi):
    return random.random()

#=======================================================================================================================
# Heleper Simulators
#=======================================================================================================================
def draw_from_prior(a, b, delta_v):
    size = len(delta_v)
    
#    alpha = 1
    alpha = gamma_rvs(a, b)
    
    labels, values = sample_from_crp(alpha, size, random.random)
    
    partition = Partition()
    
    for k, value in enumerate(values):
        partition.add_cell(value)
        
        for item, cell_index in enumerate(labels):
            if cell_index == k:
                partition.add_item(item, cell_index)

    return alpha, partition

def draw_data(partition, mu_r, mu_v, delta_v):
    g = []
    
    for d_v in delta_v:
        pi = dirichlet_rvs(d_v)
        
        g.append(discrete_rvs(pi))
     
    d = 100
    
    data = []
    
    phi = partition.item_values
    
    for i in range(len(phi)):
        mu = (1 - phi[i]) * mu_r + phi[i] * mu_v[g[i]]
        
        a = binomial_rvs(d, mu)
        
        data_point = PyCloneData(a,
                                 d,
                                 {mu_r : 1},
                                 dict(zip(mu_v, delta_v[i])))
        
        data.append(data_point)
    
    return data
    
def marginal_conditional_simulator(a, b, mu_r, mu_v, delta_v, num_iters): 
    params = []
    
    for _ in range(num_iters):
        alpha, partition = draw_from_prior(a, b, delta_v) 
        
        params.append({'alpha' : alpha, 'phi' : partition.item_values})
        
        draw_data(partition, mu_r, mu_v, delta_v)
    
    return params

def successive_conditional_simulator(a, b, mu_r, mu_v, delta_v, num_iters):
    params_sample = []
    
    alpha, partition = draw_from_prior(a, b, delta_v)
    
    for _ in range(num_iters):
        data = draw_data(partition, mu_r, mu_v, delta_v)
        
#        new_alpha = 1
        new_partition = sample_partition(alpha, data, partition, random.random)
        
        new_alpha = sample_alpha(alpha,
                                 new_partition.number_of_cells,
                                 new_partition.number_of_items,
                                 a=a,
                                 b=b)
        
#        new_partition = partition     
        
        
        sample_cell_values(data, new_partition, uniform_proposal)
        
        params = {'alpha' : new_alpha, 'phi' : new_partition.item_values}
        
        params_sample.append(params)
        
        alpha = new_alpha
        
        partition = new_partition
    
    return params_sample
    
if __name__ == "__main__":
    delta_v = [
               [1, 1],
               [1, 1]
               ]

    mu_r = 0.999    
    mu_v = [0.001, 0.5]
    
    a = 1
    b = 1
    
    burnin = 90000
    num_iters = 100000
    thin = 100
    
    params_1 = marginal_conditional_simulator(a, b, mu_r, mu_v, delta_v, num_iters)
    params_2 = successive_conditional_simulator(a, b, mu_r, mu_v, delta_v, num_iters)    
    
    trace_1 = [x['alpha'] for x in params_1][burnin::thin]
    trace_2 = [x['alpha'] for x in params_2][burnin::thin]
    
    print compare_trace(trace_1, trace_2, lambda x: x)
    print compare_trace(trace_1, trace_2, lambda x: x ** 2)
    
    for i in range(len(delta_v)):
        trace_1 = [x['phi'][i] for x in params_1][burnin::thin]
        trace_2 = [x['phi'][i] for x in params_2][burnin::thin]
        
        print compare_trace(trace_1, trace_2, lambda x: x)
        print compare_trace(trace_1, trace_2, lambda x: x ** 2)

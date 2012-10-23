'''
Created on 2012-09-15

@author: andrew
'''
from __future__ import division

from random import gammavariate as gamma_rvs

from pyclone.run import get_log_mix_weights
from pyclone.sampler import PyCloneData, DirichletProcessSampler

from pydp.base_measures import BetaBaseMeasure
from pydp.partition import Partition
from pydp.rvs import binomial_rvs, discrete_rvs, dirichlet_rvs
from pydp.tests.simulators import sample_from_crp
from pydp.tests.geweke import compare_trace

def draw_from_prior(a, b, base_measure, delta_v):
    size = len(delta_v)
    
    alpha = gamma_rvs(a, b)
    
    partition = sample_from_crp(alpha, size, base_measure)

    return alpha, partition

def draw_data(partition, mu_r, mu_v, delta_v):
    g = []
    
    for d_v in delta_v:
        pi = dirichlet_rvs(d_v)
        
        g.append(discrete_rvs(pi))
     
    d = 100
    
    data = []
    
    phi = [params.x for params in partition.item_values]
    
    for i in range(len(phi)):
        mu = (1 - phi[i]) * mu_r + phi[i] * mu_v[g[i]]
        
        a = binomial_rvs(d, mu)
        
        log_pi_r = get_log_mix_weights([1, ])
        log_pi_v = get_log_mix_weights(delta_v[i])
        
        data_point = PyCloneData(a, d, tuple([mu_r, ]), tuple(mu_v), tuple(log_pi_r), tuple(log_pi_v))
        
        data.append(data_point)
    
    return data
    
def marginal_conditional_simulator(a, b, mu_r, mu_v, delta_v, num_iters): 
    params = []
    
    base_measure = BetaBaseMeasure(1, 1)
    
    for _ in range(num_iters):
        alpha, partition = draw_from_prior(a, b, base_measure, delta_v) 
        
        params.append({
                       'alpha' : alpha,
                       'cellular_frequencies' : [value.x for value in partition.item_values]
                       })
        
        draw_data(partition, mu_r, mu_v, delta_v)
    
    return params

def successive_conditional_simulator(a, b, mu_r, mu_v, delta_v, num_iters):
    params = []
    
    base_measure = BetaBaseMeasure(1, 1)
    
    alpha, partition = draw_from_prior(a, b, base_measure, delta_v)
    
    sampler = DirichletProcessSampler()
    
    data = draw_data(partition, mu_r, mu_v, delta_v)
    
    sampler.partition = partition
    
    for _ in range(num_iters):
        sampler.interactive_sample(data)

        params.append(sampler.state)
        
        data = draw_data(partition, mu_r, mu_v, delta_v)
    
    return params
    
if __name__ == "__main__":
    delta_v = [
               [1, ],
               [1, ],
               [1, ],
               [1, ],
               [1, ],
               [1, ],
               [1, ],
               [1, ],
               [1, ],
               [1, ]
               ]

    mu_r = 0.999    
    mu_v = [0.001, ]
    
    a = 1
    b = 1
    
    num_iters = int(2e6)
    burnin = int(1e6)
    thin = 100
    
    params_1 = marginal_conditional_simulator(a, b, mu_r, mu_v, delta_v, num_iters)
    params_2 = successive_conditional_simulator(a, b, mu_r, mu_v, delta_v, num_iters)    
    
    trace_1 = [x['alpha'] for x in params_1][burnin::thin]
    trace_2 = [x['alpha'] for x in params_2][burnin::thin]
    
#    print compare_trace(trace_1, trace_2, lambda x: x)
#    print compare_trace(trace_1, trace_2, lambda x: x ** 2)
    
    for i in range(len(delta_v)):
        trace_1 = [x['cellular_frequencies'][i] for x in params_1][burnin::thin]
        trace_2 = [x['cellular_frequencies'][i] for x in params_2][burnin::thin]
        
        print compare_trace(trace_1, trace_2, lambda x: x)
        print compare_trace(trace_1, trace_2, lambda x: x ** 2)

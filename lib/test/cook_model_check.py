from __future__ import division

import csv

from collections import defaultdict

from pyclone.run import get_log_mix_weights
from pyclone.sampler import DirichletProcessSampler, PyCloneData, BetaBaseMeasure

from pydp.rvs import binomial_rvs, dirichlet_rvs, discrete_rvs, gamma_rvs
from pydp.stats import chi_square_cdf, inverse_normal_cdf
from pydp.tests.simulators import sample_from_crp

def draw_parameters_from_prior(base_measure, size):
    alpha = gamma_rvs(1, 1)
    
    partition = sample_from_crp(alpha, size, base_measure)
    
    return alpha, partition

def draw_data(partition, mu_r, delta_r, mu_v, delta_v, d=1000):
    data = []
    
    for delta, phi in zip(delta_v, partition.item_values):
        pi = dirichlet_rvs(delta)
        
        g = discrete_rvs(pi)
        
        mu = (1 - phi.x) * mu_r[0] + phi.x * mu_v[g]
        
        a = binomial_rvs(d, mu)
        
        log_pi_r = get_log_mix_weights(delta_r)
        log_pi_v = get_log_mix_weights(delta)
        
        data.append(PyCloneData(a, d, mu_r, mu_v, log_pi_r, log_pi_v))
    
    return data

def compute_quartiles(alpha, partition, sampler_states):
    quartiles = {}
    
    trace = [x['alpha'] for x in sampler_states]
    
    quartiles['alpha'] = _compute_quartile(alpha, trace)
    
    for i, true_phi in enumerate(partition.item_values):
        trace = [x['cellular_frequencies'][i] for x in sampler_states]
        
        quartiles['p_{0}'.format(i)] = _compute_quartile(true_phi.x, trace)
    
    return quartiles

def _compute_quartile(true_value, trace):
    q = 0
    
    for value in trace:
        if true_value > value:
            q += 1
    
    q = q / len(trace)
    
#    print true_value, mean(trace), q, inverse_normal_cdf(q) 
    
    return q

#def compute_significance(test_statistics):
#    for param in test_statistics:
#        test_stat = test_statistics[param]
#        
#        inverse_normal_cdf(q) ** 2
#        
#        print param, 1 - chi_square_cdf(sum(test_stat), len(test_stat))
    
#=======================================================================================================================
# Simulators
#=======================================================================================================================
size = 4

num_iters = int(1e5)
burnin = int(9e4)
thin = 1
num_replicates = 100

# Pyclone Parameters
eps = 1e-3

mu_r = [1 - eps, ]

delta_r = [1, ]

mu_v = [eps, 0.5]

delta_v = [[1, 1]] * size

# Run Test
base_measure = BetaBaseMeasure(1, 1)

quartiles = defaultdict(list)

for i in range(num_replicates):
    print i
    
    alpha, partition = draw_parameters_from_prior(base_measure, size)
    
    observed_p = alpha
    
    params = [param.x for param in partition.item_values]
    
    data = draw_data(partition, mu_r, delta_r, mu_v, delta_v)
    
    sampler = DirichletProcessSampler()
    
    sampler.initialise_partition(data)
    
    sampler_states = []
    
    for _ in range(num_iters):
        sampler.interactive_sample(data)

        sampler_states.append(sampler.state)
    
    replicate_quartiles = compute_quartiles(alpha, partition, sampler_states)
    
    for param in replicate_quartiles:
        quartiles[param].append(replicate_quartiles[param])

writer = csv.writer(open('cook_test.tsv', 'w'), delimiter='\t')

for param in quartiles:
    out_row = [param, ] + quartiles[param]

    writer.writerow(out_row)
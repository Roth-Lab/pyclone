'''
Created on 2012-09-15

@author: andrew
'''
from __future__ import division

import numpy as np
import random

from random import gammavariate as gamma_rvs

from pyclone.utils import discrete_rvs

from pyclone._likelihoods import BinomialLikelihood
#from pyclone._sampler import Customer, Restaurant

from pyclone.test.crp import ChineseRestaurantProcess
from pyclone.test.geweke import compare_trace


#=======================================================================================================================
# Heleper Simulators
#=======================================================================================================================
def draw_from_prior(priors):
    params = {}
    
    n = priors['delta_v'].shape[0]

    #params['alpha'] = 1
    params['alpha'] = gamma_rvs(priors['alpha']['location'], priors['alpha']['scale'])
    
    crp = ChineseRestaurantProcess(params['alpha'], n, random.random)
    
    params['phi'] = crp.draw_dishes()
    
    params['mu_r'] = priors['mu_r']
    params['mu_v'] = priors['mu_v']
    
    params['delta_v'] = priors['delta_v']
    
    params['d'] = np.zeros(n, dtype='int')
    params['d'][:] = 1000
    
    return params

def draw_data(params):
    g = []
        
    for x in params['delta_v']:
        pi = np.random.dirichlet(x)
        
        g.append(discrete_rvs(pi))
   
    mu_v = params['mu_v'][g]
    
    mu_r = params['mu_r']
    
    phi = params['phi']
    
    mu = (1 - phi) * mu_r + phi * mu_v
    
    a = np.random.binomial(params['d'], mu)
    
    return a    
    
def marginal_conditional_simulator(priors, num_iters):
    data = []    
    params = []
    
    for _ in range(num_iters):
        params.append(draw_from_prior(priors))
    
        data.append(draw_data(params[-1]))
    
    return data, params

def successive_conditional_simulator(priors, num_iters):
    data = []
    params = []
    
    restaurant = None
    
    old_params = draw_from_prior(priors)  
    
    for _ in range(num_iters):
        a = draw_data(old_params)
        
        d = old_params['d']
        
        mu_r = priors['mu_r']
        ref_priors = {mu_r : 1}
                    
        likelihoods = []
        
        for i, (a_i, d_i) in enumerate(zip(a, d)):
            var_priors = {}
            
            for mu_v, delta_v in zip(priors['mu_v'], priors['delta_v'][i]):
                var_priors[mu_v] = delta_v                           
            
            likelihoods.append(BinomialLikelihood(a_i, d_i, ref_priors, var_priors))
        
        if restaurant is None:
            customers = [Customer(x) for x in likelihoods]
            
            restaurant = Restaurant(customers, concentration=old_params['alpha'], m=2)
        else:
            for customer, likelihood in zip(restaurant.customers, likelihoods):
                customer._likelihood = likelihood

        restaurant.sample_partition()
        
        restaurant.sample_cellular_frequencies()
        
        restaurant.sample_concentration()
        
        new_params = old_params.copy()
        
        new_params['alpha'] = restaurant.concentration
        
        phi = []
        
        for table_id in restaurant.seating_assignment:
            phi.append(restaurant.dishes[table_id])
        
        new_params['phi'] = np.array(phi)
        
        params.append(new_params)            
        data.append(a)
        
        old_params = new_params
    
    return data, params

if __name__ == "__main__":
    priors = {}
    
    priors['delta_v'] = np.array([
                                  [1,],
                                  [1,]   
                                  ])

    priors['mu_r'] = 0.999    
    priors['mu_v'] = np.array([0.5,])
    
    priors['depth'] = 1000
    
    priors['alpha'] = {}
    priors['alpha']['location'] = 1
    priors['alpha']['scale'] = 1
    
    burnin = 90000
    num_iters = 100000
    thin = 1
    
    data_1, params_1 = marginal_conditional_simulator(priors, num_iters)
    data_2, params_2 = marginal_conditional_simulator(priors, num_iters)    
    
    trace_1 = np.array([x['alpha'] for x in params_1][burnin::thin])
    trace_2 = np.array([x['alpha'] for x in params_2][burnin::thin])
    
    print compare_trace(trace_1, trace_2, lambda x: x)
    print compare_trace(trace_1, trace_2, lambda x: x**2)
    
    for i in range(len(priors['delta_v'])):
        trace_1 =  np.array([x['phi'][i] for x in params_1][burnin::thin])
        trace_2 =  np.array([x['phi'][i] for x in params_2][burnin::thin])
        
        print compare_trace(trace_1, trace_2, lambda x: x)
        print compare_trace(trace_1, trace_2, lambda x: x**2)

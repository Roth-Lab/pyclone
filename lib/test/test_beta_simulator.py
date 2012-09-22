'''
Created on 2012-09-15

@author: andrew
'''
import numpy as np
import random

from random import gammavariate as gamma_rvs

from pyclone.utils import discrete_rvs

from pyclone._likelihoods import BetaBinomialLikelihood
from pyclone._sampler import Customer, Restaurant


#=======================================================================================================================
# Heleper Simulators
#=======================================================================================================================
def draw_from_prior(priors):
    params = {}
    
    n = priors['delta_v'].shape[0]

    params['alpha'] = gamma_rvs(priors['alpha']['location'], priors['alpha']['scale'])
    
    crp = ChineseRestaurantProcess(params['alpha'], n)
    
    params['phi'] = crp.draw_phi()
    
    params['mu_r'] = priors['mu_r']
    
    params['delta_v'] = priors['delta_v']
    
    params['d'] = np.random.poisson(priors['depth'], n)
    
    s = 100#random.uniform(1, 10000)
    mu_v = []
    
    for m in priors['mu_v']:
        a = s * m
        b = s - a
        
        mu_v.append(random.betavariate(a, b))
    
    params['s'] = s
    params['mu_v'] = np.array(mu_v) 
    
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
            
            for m, d in zip(priors['mu_v'], priors['delta_v'][i]):
                var_priors[m] = d                                
            
            likelihoods.append(BetaBinomialLikelihood(a_i, d_i, ref_priors, var_priors))
        
        if restaurant is None:
            customers = [Customer(x) for x in likelihoods]
            
            restaurant = Restaurant(customers, concentration=old_params['alpha'], beta_precision=100)
        else:
            for customer, likelihood in zip(restaurant.customers, likelihoods):
                customer._likelihood = likelihood
        
        restaurant.sample_partition()
        
        restaurant.sample_cellular_frequencies()
        
        restaurant.sample_concentration()
        
#        restaurant.sample_beta_precision()
        
        new_params = old_params.copy()
        
        new_params['alpha'] = restaurant.concentration
        
        phi = []
        
        for table_id in restaurant.seating_assignment:
            phi.append(restaurant.dishes[table_id])
        
        new_params['phi'] = np.array(phi)
        
        new_params['s'] = restaurant.beta_precision
        
        params.append(new_params)            
        data.append(a)
        
        old_params = new_params
        
    
    return data, params

class ChineseRestaurantProcess(object):
    def __init__(self, alpha, num_customers):
        self.alpha = alpha
        
        self.num_customers = num_customers
        
        self._seat_customers()
    
    def draw_phi(self):
        phi = []
        
        for i in range(self.num_customers):
            for table, dish in zip(self.tables, self.values):            
                if i in table:
                    phi.append(dish)
                
        return np.array(phi)
    
    def _seat_customers(self):
        tables = []
        values = []
        
        # Seat the first customer
        tables.append([0, ])
        values.append(random.random())
        
        for customer in range(1, self.num_customers):
            p = self._get_table_probabilities(tables, self.alpha)
            
            table_id = discrete_rvs(p)
            
            if table_id == len(tables):
                tables.append([customer, ])
                values.append(random.random())
            else:
                tables[table_id].append(customer)
        
        self.tables = tables
        self.values = values
                        
    def _get_table_probabilities(self, tables, alpha):
        p = []
        
        for table in tables:
            p.append(len(table))
        
        p.append(alpha)
        
        p = [x / sum(p) for x in p]
        
        return p

if __name__ == "__main__":
    priors = {}
    
    priors['delta_v'] = np.array([
                                  [1, 1],
                                  [1, 1]   
                                  ])

    priors['mu_r'] = 0.999    
    priors['mu_v'] = np.array([0.001, 0.5])
    
    priors['depth'] = 1000
    
    priors['alpha'] = {}
    priors['alpha']['location'] = 1
    priors['alpha']['scale'] = 1
    
    burnin = 1
    num_iters = 100
    
    data_1, params_1 = marginal_conditional_simulator(priors, num_iters)
    data_2, params_2 = successive_conditional_simulator(priors, num_iters)

    print sum([x['phi'][0] for x in params_1][burnin::]) / (num_iters - burnin)
    print sum([x['phi'][0] for x in params_2][burnin::]) / (num_iters - burnin)

    print sum([x['phi'][1] for x in params_1][burnin::]) / (num_iters - burnin)
    print sum([x['phi'][1] for x in params_2][burnin::]) / (num_iters - burnin)

    print sum([x['phi'][0]**2 for x in params_1][burnin::]) / (num_iters - burnin)
    print sum([x['phi'][0]**2 for x in params_2][burnin::]) / (num_iters - burnin)

    print sum([x['phi'][1]**2 for x in params_1][burnin::]) / (num_iters - burnin)
    print sum([x['phi'][1]**2 for x in params_2][burnin::]) / (num_iters - burnin)
    
    print sum([x['alpha'] for x in params_1][burnin::]) / (num_iters - burnin)
    print sum([x['alpha'] for x in params_2][burnin::]) / (num_iters - burnin)
    
    print sum([x['alpha'] ** 2 for x in params_1][burnin::]) / (num_iters - burnin)
    print sum([x['alpha'] ** 2 for x in params_2][burnin::]) / (num_iters - burnin)
    
    print sum(data_1[burnin::]) / (num_iters - burnin)
    print sum(data_2[burnin::]) / (num_iters - burnin)

    print sum([x['s'] for x in params_1][burnin::]) / (num_iters - burnin)
    print sum([x['s'] for x in params_2][burnin::]) / (num_iters - burnin)
'''
Created on 2012-06-14

@author: Andrew Roth
'''
from __future__ import division

import random

from numpy.random import dirichlet
from pyclone.utils import discrete_rvs

class ChineseRestaurantProcess(object):
    def __init__(self, alpha, num_customers):
        self.alpha = alpha
        
        self.num_customers = num_customers
        
        self._current_customer = 0
        
        self._seat_customers()
    
    def draw_cellular_frequency(self):
        if self._current_customer >= self.num_customers:
            raise Exception("To many draws from CRP.")
        
        for table, phi in zip(self.tables, self.values):
            if self._current_customer in table:
                self._current_customer += 1
                
                return phi
    
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

class MixtureModel(object):
    '''
    Object for simulating cellular frequencies.        
    '''    
    def __init__(self, num_clones, mix_weight_priors=None, frequencies=None):
        '''
        Args:
            num_clones (int): Number of clonal populations.
        
        Kwargs:
            mix_weight_priors (list ints): Dirichlet prior pseudo counts for each class.
            frequencies (list floats) : Frequencies for each class to be assosciated with each class.
        ''' 
        self.num_clones = num_clones
        
        if mix_weight_priors is None:
            mix_weight_priors = [2., ] * self.num_clones
        else:
            if len(mix_weight_priors) != num_clones:
                raise Exception("Length of mix-weight priors must equal number of clonal populations.")
        
        self.mix_weights = dirichlet(mix_weight_priors)
        
        if frequencies is None:
            self.frequencies = [random.random() for _ in range(num_clones)]
        else:
            if len(frequencies) != num_clones:
                raise Exception("Number of frequencies passed in must match the number of clonal populations.")
            
            self.frequencies = frequencies
    
    def draw_cellular_frequency(self):
        '''
        Draw a cellular frequency.
        '''
        population_id = discrete_rvs(self.mix_weights)
        
        return self.frequencies[population_id]
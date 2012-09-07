'''
Created on 2011-12-15

@author: Andrew Roth
'''
from __future__ import division

from math import exp, log

from random import betavariate as beta_rvs, gammavariate as gamma_rvs, normalvariate as normal_rvs, \
    uniform as uniform_rvs

from pyclone.utils import discrete_rvs, log_space_normalise

class DirichletProcessSampler(object):
    def __init__(self, likelihoods, model, concentration=None, m=2):
        customers = [Customer(x) for x in likelihoods]
        
        if concentration is None:
            self.update_concentration = True
            
            concentration = 1
        else:
            self.update_concentration = False
        
        if model == 'binomial':
            self.update_beta_precision = False            
            beta_precision = None
            
        elif model == 'beta-binomial':
            self.update_beta_precision = True            
            beta_precision = 100
                        
        else:
            raise Exception("Model {0} not recognised.".format(model))

        self.restaurant = Restaurant(customers, beta_precision=beta_precision, concentration=concentration, m=m)
        
        self.num_iters = 0
        
    def sample(self, results_db, num_iters, print_freq=100):
        for _ in range(num_iters):
            if self.num_iters % print_freq == 0:
                print self.num_iters, self.restaurant.number_of_tables, self.restaurant.concentration,
                
                if self.update_beta_precision:
                    print self.restaurant.beta_precision
                else:
                    print
            
            self.num_iters += 1
            
            self.restaurant.sample_partition()
            
            self.restaurant.sample_cellular_frequencies()
            
            if self.update_concentration:
                self.restaurant.sample_concentration()
            
            if self.update_beta_precision:
                self.restaurant.sample_beta_precision()
            
            results_db.update_trace(self.state)
                
    @property
    def state(self):
        state = {
                  'alpha' : self.restaurant.concentration,
                  'labels' : self.restaurant.seating_assignment,
                  'phi' : self.restaurant.dishes                
                  }
        
        if self.update_beta_precision:
            state['beta_precision'] = self.restaurant.beta_precision
        
        return state

class Restaurant(object):
    def __init__(self, customers, beta_precision=None, concentration=1, m=2):
        self.customers = customers
        
        self.beta_precision = beta_precision
        
        self.concentration = concentration
        
        self.m = m
        
        self._init_tables()

    @property
    def dishes(self):
        '''
        Returns a list of dishes served at each table
        '''
        dishes = []
        
        for table in self.tables:
            dishes.append(table.dish[0])
        
        return dishes
    
    @property
    def number_of_tables(self):
        return len(self.tables)

    @property
    def seating_assignment(self):
        '''
        Returns a list with the id of the table each customer is seated at.
        '''        
        seating_assignment = []
        
        for customer in self.customers:
            table = self.table_map[customer]
            
            table_id = self.tables.index(table)
            
            seating_assignment.append(table_id)
        
        return seating_assignment
                
    def sample_beta_precision(self):
        '''
        Metropolis-Hastings update of beta precision.
        '''
        old_beta_precision = self.beta_precision
        
        eps = normal_rvs(0, 1)
        
        new_beta_precision = exp(log(old_beta_precision) + eps)

        old_ll = 0
        new_ll = 0
        
        for table in self.tables:
            phi = table.dish[0]
            
            old_dish = phi, old_beta_precision
            new_dish = phi, new_beta_precision
            
            old_ll += table.evaluate_table_likelihood(old_dish)
            new_ll += table.evaluate_table_likelihood(new_dish)
        
        log_ratio = new_ll - old_ll
        
        u = uniform_rvs(0, 1)
        
        if log_ratio >= log(u):
            self.beta_precision = new_beta_precision
        
            for table in self.tables:
                phi = table.dish[0]
                
                table.dish = phi, new_beta_precision
   
    def sample_cellular_frequencies(self):
        '''
        Metropolis-Hastings update of cellular frequencies for each table.
        '''
        for table in self.tables:
            old_phi = table.dish[0]

            if old_phi == 1:
                old_phi = 1 - 1e-10
            elif old_phi == 0:
                old_phi = 1e-10 
    
            old_phi_star = old_phi / (1 - old_phi)
            
            eps = normal_rvs(0, 1)
            
            new_phi_star = exp(log(old_phi_star) + eps)
            
            new_phi = new_phi_star / (1 + new_phi_star) 
            
            if self.beta_precision is None:
                new_dish = (new_phi,) 
            else:
                new_dish = (new_phi, self.beta_precision)
            
            table.try_new_dish(new_dish)

    def sample_concentration(self):
        a = 1
        b = 1
        
        n = len(self.customers)
        k = len(self.tables)
        
        eta = beta_rvs(self.concentration + 1, n)

        x = (a + k - 1) / (n * (b - log(eta)))
        
        pi = x / (1 + x)

        label = discrete_rvs([pi, 1 - pi])
        
        scale = 1 / (b - log(eta))
                
        if label == 0:
            self.concentration = gamma_rvs(a + k, scale)
        elif label == 1:
            self.concentration = gamma_rvs(a + k - 1, scale)

    def sample_partition(self):
        for customer in self.customers:
            self._reseat_customer(customer)
    
    def _reseat_customer(self, customer):
        '''
        Reseat customer using algorithm 8 from 
        Markov chain sampling methods for Dirichlet process mixture models, R.M. Neal, 2000.
        '''
        table = self.table_map[customer]
        
        table.remove_customer(customer)
        
        # Create auxillary tables
        if table.empty:
            num_new_tables = self.m - 1
        else:
            num_new_tables = self.m 
        
        for _ in range(num_new_tables):
            table = self._create_table()
            
            self.tables.append(table)
        
        # Sample new table for the customer
        log_p = []
        
        for table in self.tables:
            dish = table.dish
            
            ll = customer.likelihood.evaluate(*dish)
            
            if table.empty:
                counts = self.concentration / self.m
            else:
                counts = table.number_of_customers
            
            log_p.append(log(counts) + ll)
            
        log_p = log_space_normalise(log_p)
        
        p = [exp(x) for x in log_p]
        
        table_id = discrete_rvs(p)
        
        # Assign customer to table update map
        table = self.tables[table_id]
        
        table.add_customer(customer)
        
        self.table_map[customer] = table
        
        # Remove empty tables
        empty_tables = []
        
        for table in self.tables:
            if table.empty:
                empty_tables.append(table)
        
        for table in empty_tables:
            self.tables.remove(table)

    def _init_tables(self):
        self.tables = []
        self.table_map = {}
        
        for customer in self.customers:
            table = self._create_table()
            
            table.add_customer(customer)
            
            self.tables.append(table)
            
            self.table_map[customer] = table
        
        # Initialise clusters to good values by re-sampling
        for _ in range(100):
            self.sample_cellular_frequencies()
    
    def _create_table(self):
        table = Table()

        new_phi = self._cellular_frequency_proposal() 
        
        if self.beta_precision is None:
            new_dish = (new_phi,) 
        else:
            new_dish = (new_phi, self.beta_precision)
            
        table.dish = new_dish
        
        return table

    def _cellular_frequency_proposal(self):
        return uniform_rvs(0, 1)

class Table(object):
    def __init__(self):
        self._dish = None
        self._customers = []
    
    @property
    def customers(self):
        return self._customers
    
    @property
    def dish(self):
        return self._dish
    
    @dish.setter
    def dish(self, value):
        self._dish = value

    @property
    def empty(self):
        if len(self._customers) == 0:
            return True
        else:
            return False
    
    @property
    def number_of_customers(self):
        return len(self._customers)
        
    def add_customer(self, customer):
        self._customers.append(customer)
    
    def remove_customer(self, customer):
        self._customers.remove(customer)
    
    def try_new_dish(self, new_dish):
        old_dish = self.dish
        
        old_ll = self.evaluate_table_likelihood(old_dish)
        new_ll = self.evaluate_table_likelihood(new_dish)
        
        log_ratio = new_ll - old_ll
        
        u = uniform_rvs(0, 1)
        
        if log_ratio >= log(u):
            self.dish = new_dish
            
    def evaluate_table_likelihood(self, dish):
        ll = 0
        
        for customer in self._customers:
            ll += customer.likelihood.evaluate(*dish)
        
        return ll

class Customer(object):
    def __init__(self, likelihood):
        self._likelihood = likelihood
    
    @property   
    def likelihood(self):
        return self._likelihood

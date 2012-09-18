'''
Created on 2011-12-15

@author: Andrew Roth
'''
from __future__ import division

from math import exp, log

from random import betavariate as beta_rvs, gammavariate as gamma_rvs, uniform as uniform_rvs, shuffle

from pyclone.utils import discrete_rvs, log_space_normalise

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

def sample_phi(old_value, data, priors):
        
class DirichletProcessSampler(object):
    def __init__(self, customers, concentration_parameter, update_concentration=True, m=2):
        self.customers = customers
        
        self.cook = UniformCook(0, 1)
        
        self.host = Host(m)
        
        self.restaurant = Restaurant(self.cook, customers)
        
        self.concentration_parameter = concentration_parameter
        
        self.update_concentration = update_concentration
        
        self.num_iters = 0
        
    def sample(self, results_db, num_iters, print_freq=100):
        for _ in range(num_iters):
            if self.num_iters % print_freq == 0:
                print self.num_iters, self.restaurant.number_of_tables, self.concentration_parameter.value
            
            state = self.interactive_sample()
            
            results_db.update_trace(state)
    
    def interactive_sample(self):
        self.num_iters += 1
        
        if self.update_concentration:
            self._sample_concentration_parameter()
        
        self._sample_partition()
        
        self._sample_cluster_parameters()
        
        return self.state
    
    def _sample_concentration_parameter(self):
        self.concentration_parameter.sample(self.restaurant.number_of_tables,
                                            self.restaurant.number_of_customers)
    
    def _sample_partition(self):
        self.host.reseat_customers(self.concentration.value,
                                   self.restaurant)
    
    def _sample_cluster_parameters(self):
        self.restaurant.try_new_dishes()
                
    @property
    def state(self):
        state = {
                 'alpha' : self.concentration_parameter.value,
                 'labels' : self.restaurant.seating_assignment,
                 'phi' : self.restaurant.dishes                
                 }

        return state

class ConcentrationParameter(object):
    def __init__(self, value, a=1, b=1):
        '''
        Args:
            value : (float) Initial value.
        
        Kwargs:
            a : (float) Location parameter in gamma prior
            b : (float) Scale parameter in gamma prior
        '''
        self.value = value
        
        self.a = a
        self.b = b
        
    def sample(self, k, n):
        '''
        Args:
            k : (int) Number of active clusters in DP (number of tables in restaurant)
            n : (int) Number of data points in DP (number of customers in restaurant)
        '''
        a = self.a
        b = self.b
        
        old_value = self.value
        
        eta = beta_rvs(old_value + 1, n)

        x = (a + k - 1) / (n * (b - log(eta)))
        
        pi = x / (1 + x)

        label = discrete_rvs([pi, 1 - pi])
        
        scale = 1 / (b - log(eta))
                
        if label == 0:
            new_value = gamma_rvs(a + k, scale)
        else:
            new_value = gamma_rvs(a + k - 1, scale)
        
        self.value = new_value
        
class Host(object):
    '''
    Class responsible for reseating customers in restaurant.
    '''
    def __init__(self, m=2):
        self.m = m
    
    def reseat_customers(self, alpha, restaurant):
        customers = restaurant.customers
        
        shuffle(customers)
        
        for customer in customers:
            self._reseat_customer(alpha, customer, restaurant)
        
    def _reseat_customer(self, alpha, customer, restaurant):
        '''
        Reseat customer using algorithm 8 from Markov chain sampling methods for Dirichlet process mixture models, R.M.
        Neal, 2000.
        '''
        tables = restaurant.tables
        
        table_map = restaurant.table_map
        
        table = table_map[customer]
        
        table.remove_customer(customer)
        
        # Create auxillary tables
        if table.empty:
            num_new_tables = self.m - 1
        else:
            num_new_tables = self.m
        
        for _ in range(num_new_tables):
            tables.append(Table(self.cook))
        
        # Sample new table for the customer
        log_probs = []
        
        for table in tables:
            cluster_log_prob = customer.logp(table.dish)
            
            if table.empty:
                counts = alpha / self.m
            else:
                counts = table.number_of_customers
            
            log_probs.append(log(counts) + cluster_log_prob)
            
        log_probs = log_space_normalise(log_probs)
        
        p = [exp(x) for x in log_probs]
        
        table_id = discrete_rvs(p)
        
        # Assign customer to table update map
        table = tables[table_id]
        
        table.add_customer(customer)
        
        table_map[customer] = table
        
        # Remove empty tables
        for table in tables[:]:
            if table.empty:
                tables.remove(table)
        
        restaurant.tables = tables
        restaurant.table_map = table_map    
            
        
class Restaurant(object):
    def __init__(self, cook, customers):
        self.cook = cook
        
        self._customers = customers
        
        self._init_tables()

    @property
    def customers(self):
        '''
        Return a list of customers in restaurant.
        '''
        return self._customers[:]

    @property
    def dishes(self):
        '''
        Returns a list of values assosciated with each cluster (dishes served at each table).
        '''
        dishes = []
        
        for table in self._tables:
            dishes.append(table.dish[0])
        
        return dishes
    
    @property
    def number_of_customers(self):
        '''
        Returns the number of data points (customers) currently used by DP.
        '''
        return len(self._customers)
    
    @property
    def number_of_tables(self):
        '''
        Returns the number of clusters (tables) currently used by the DP.
        '''
        return len(self._tables)

    @property
    def seating_assignment(self):
        '''
        Returns a list with the id of the table each customer is seated at.
        '''        
        seating_assignment = []
        
        for customer in self._customers:
            table = self._table_map[customer]
            
            table_id = self._tables.index(table)
            
            seating_assignment.append(table_id)
        
        return seating_assignment
    
    @property
    def tables(self):
        return self._tables[:]
    
    @tables.setter
    def tables(self, value):
        self._tables = value[:]
        
    @property
    def table_map(self):
        return self._table_map.copy()
    
    @table_map.setter
    def table_map(self, value):
        self._table_map = value.copy()
   
    def try_new_dishes(self):
        '''
        Metropolis-Hastings update of cellular frequencies for each table.
        '''
        for table in self.tables:
            table.try_new_dish()

    def _init_tables(self):
        self._tables = []
        self._table_map = {}
        
        for customer in self.customers:
            table = Table(self.cook)
            
            table.add_customer(customer)
            
            self._tables.append(table)
            
            self._table_map[customer] = table

class Table(object):
    def __init__(self, cook):
        self.cook = cook
        
        self._dish = cook.propose_dish()
        
        self._customers = []
    
    @property
    def customers(self):
        return self._customers[:]
    
    @property
    def dish(self):
        return self._dish
    
    @dish.setter
    def dish(self, value):
        self._dish = value

    @property
    def empty(self):
        if self.number_of_customers == 0:
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
    
    def try_new_dish(self):
        new_dish = self.cook.propose_dish()
        
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
            ll += customer.logp(dish)
        
        return ll

class UniformCook(object):
    '''
    Cook for proposing dishes from Uniform(a,b) continuous distribution.
    '''
    def __init__(self, a=0, b=1):
        '''
        Kwargs:
            a : (float) Left end point of Uniform distribution
            b : (float) Right end point of Uniform distribution
        ''' 
        self.a = a
        self.b = b
    
    def propose_dish(self):
        return uniform_rvs(self.a, self.b)

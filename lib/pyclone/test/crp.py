'''
Created on 2012-09-18

@author: andrew
'''
from pyclone.utils import discrete_rvs

import numpy as np

class ChineseRestaurantProcess(object):
    def __init__(self, alpha, num_customers, proposal_func):
        self.alpha = alpha
        
        self.num_customers = num_customers
        
        self.proposal_func = proposal_func
        
        self._seat_customers()
    
    def draw_dishes(self):
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
        values.append(self.proposal_func())
        
        for customer in range(1, self.num_customers):
            p = self._get_table_probabilities(tables, self.alpha)
            
            table_id = discrete_rvs(p)
            
            if table_id == len(tables):
                tables.append([customer, ])
                values.append(self.proposal_func())
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

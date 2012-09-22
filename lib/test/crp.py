'''
Created on 2012-09-18

@author: andrew
'''
from __future__ import division

from pyclone.utils import discrete_rvs

def sample_from_crp(alpha, size, base_measure):
    labels = []
    values = []
    
    tables = []
    
    # Seat the first customer
    tables.append([0, ])
    
    labels.append(0)
    values.append(base_measure.random())
    
    for customer in range(1, size):
        p = _get_table_probabilities(tables, alpha)
        
        table_id = discrete_rvs(p)
        
        if table_id == len(tables):
            tables.append([customer, ])
            
            values.append(base_measure.random())
        else:
            tables[table_id].append(customer)
        
        labels.append(table_id)
    
    return labels, values
                        
def _get_table_probabilities(tables, alpha):
    p = []
    
    for table in tables:
        p.append(len(table))
    
    p.append(alpha)
    
    p = [x / sum(p) for x in p]
    
    return p

'''
Created on 2012-05-15

@author: Andrew Roth
'''
from collections import OrderedDict

import csv

class DataSet(object):
    def __init__(self, file_name):
        self._data = OrderedDict()
    
        reader = csv.DictReader(open(file_name), delimiter='\t')
    
        for row in reader:
            gene = row['gene']
            
            a = int(row['a'])
            
            d = int(row['d'])
            
            mu_r = [float(x) for x in row['mu_r'].split(',')]
            delta_r = [float(x) for x in row['delta_r'].split(',')]
            
            mu_v = [float(x) for x in row['mu_v'].split(',')]
            delta_v = [float(x) for x in row['delta_v'].split(',')]
            
            self._data[gene] = DataPoint(a, d, mu_r, mu_v, delta_r, delta_v)
    
    def __iter__(self):
        return self._data.items()

class DataPoint(object):
    def __init__(self, a, d, mu_r, mu_v, delta_r, delta_v):
        self.a = a
        self.d = d
               
        self.mu_r = mu_r
        self.mu_v = mu_v
        
        self.delta_r = delta_r
        self.delta_v = delta_v 
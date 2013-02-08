'''
Utility functions for post-processing the results of a PyClone analysis.

Created on 2013-02-08

@author: Andrew Roth
'''
from collections import defaultdict

import bz2
import csv

def load_cellular_frequencies_trace(file_name, burnin, thin):
    return _load_trace(file_name, burnin, thin, float)

def load_cluster_labels_trace(file_name, burnin, thin):
    return _load_trace(file_name, burnin, thin, int)

def _load_trace(file_name, burnin, thin, cast_func):
    labels = defaultdict(list)
    
    reader = csv.DictReader(bz2.BZ2File(file_name), delimiter='\t')
    
    for i, row in enumerate(reader):
        if i < burnin:
            continue
        
        if i % thin != 0:
            continue
        
        for mutation in row:
            labels[mutation].append(cast_func(row[mutation])) 
    
    return labels

'''
Functions for clustering results from PyClone analysis.

Created on 2013-02-08

@author: Andrew Roth
'''
from __future__ import division

from pydp.cluster import cluster_with_mpear

import csv

from pyclone.post_process.utils import load_cluster_labels_trace

def write_pyclone_cluster_file(labels_file, cluster_file, burnin, thin):
    labels = cluster_pyclone_trace(labels_file, burnin, thin)
    
    writer = csv.DictWriter(
                            open(cluster_file, 'w'),
                            ['mutation_id', 'cluster_id'],
                            delimiter='\t'
                            )
    
    writer.writeheader()
    
    for mutation_id, cluster_id in labels.items():
        out_row = {'mutation_id' : mutation_id, 'cluster_id' : int(cluster_id)}
        
        writer.writerow(out_row)  
        
def cluster_pyclone_trace(labels_file, burnin, thin):    
    trace = load_cluster_labels_trace(labels_file, burnin, thin)
    
    print trace.shape
    
    labels = cluster_with_mpear(trace)
    
    mutation_ids = trace.columns
    
    return dict(zip(mutation_ids, labels))
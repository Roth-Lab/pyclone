'''
Functions for clustering results from PyClone analysis.

Created on 2013-02-08

@author: Andrew Roth
'''
from pydp.cluster import cluster_with_mpear

import pandas as pd

from pyclone.post_process.utils import load_cluster_labels_trace

def write_pyclone_cluster_file(labels_file, cluster_file, burnin, thin):
    labels = cluster_pyclone_trace(labels_file, burnin, thin)
    
    labels.to_csv(cluster_file, sep='\t')

def cluster_pyclone_trace(labels_file, burnin, thin):    
    trace = load_cluster_labels_trace(labels_file, burnin, thin)
    
    X = trace.values

    labels = cluster_with_mpear(X)
    
    labels = pd.Series(labels, index=trace.columns)
    
    labels = labels.reset_index()
    
    labels.columns = 'mutation_id', 'cluster_id'
    
    return labels

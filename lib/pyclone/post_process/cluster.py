'''
Functions for clustering results from PyClone analysis.

Created on 2013-02-08

@author: Andrew Roth
'''
from __future__ import division

import numpy as np

from math import exp
from pydp.densities import log_binomial_coefficient

try:
    import pandas as pd
except:
    raise Exception("The multi sample plotting module requires the pandas package. See http://http://pandas.pydata.org.")

try:
    from scipy.cluster.hierarchy import average, fcluster 
    from scipy.spatial.distance import pdist, squareform
except:
    raise Exception("The clustering module requires the scipy package. See http://www.scipy.org/.")

from pyclone.post_process.utils import load_trace

def write_pyclone_cluster_file(labels_file, cluster_file, burnin, thin):
    labels = cluster_pyclone_trace(labels_file, burnin, thin)
    
    labels.to_csv(cluster_file, index_label='mutation_id', sep='\t')
        
def cluster_pyclone_trace(labels_file, burnin, thin):    
    trace = load_trace(labels_file, burnin, thin)
   
    labels = cluster_with_mpear(trace)

    labels = pd.Series(labels, index=trace.index)
    
    return pd.DataFrame(labels, columns=['cluster_id'])

def cluster_with_mpear(X):
    X = np.array(X)
    
    dist_mat = pdist(X, metric='hamming')
    
    sim_mat = 1 - squareform(dist_mat)
    
    Z = average(dist_mat)
    
    max_pear = 0

    best_cluster_labels = _get_flat_clustering(Z, 1)
    
    for i in range(1, len(X) + 1):
        cluster_labels = _get_flat_clustering(Z, i)
    
        pear = _compute_mpear(cluster_labels, sim_mat)
        
        if pear > max_pear:
            max_pear = pear

            best_cluster_labels = cluster_labels
    
    return best_cluster_labels

def _get_flat_clustering(Z, number_of_clusters):
    N = len(Z) + 1
    
    if number_of_clusters == N:
        return np.arange(1, N + 1)
    
    return fcluster(Z, number_of_clusters, criterion='maxclust')

def _compute_mpear(cluster_labels, sim_mat):
    N = sim_mat.shape[0]
    
    c = exp(log_binomial_coefficient(N, 2))
    
    num_term_1 = 0
    
    for j in range(N):
        for i in range(j):
            if cluster_labels[i] == cluster_labels[j]:
                num_term_1 += sim_mat[i][j]

    num_term_2 = 0
    
    for j in range(N):
        for i in range(j):
            if cluster_labels[i] == cluster_labels[j]:
                num_term_2 += sim_mat[:j - 1, j].sum()
    
    num_term_2 /= c
    
    den_term_1 = 0
    
    for j in range(N):
        for i in range(j):
            den_term_1 += sim_mat[i][j]
            
            if cluster_labels[i] == cluster_labels[j]:
                den_term_1 += 1
    
    den_term_1 /= 2
    
    num = num_term_1 - num_term_2
    den = den_term_1 - num_term_2
    
    return num / den
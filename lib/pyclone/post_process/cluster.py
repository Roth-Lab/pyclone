'''
Functions for clustering results from PyClone analysis.

Created on 2013-02-08

@author: Andrew Roth
'''
from __future__ import division

from pydp.cluster import cluster_with_mpear

import csv
import numpy as np

try:
    from scipy.cluster.hierarchy import fclusterdata 
    from scipy.spatial.distance import pdist, squareform
    
except:
    raise Exception("The clustering module requires the scipy package. See http://www.scipy.org/.")

from pyclone.post_process.utils import load_cluster_labels_trace

def write_pyclone_cluster_file(labels_file, cluster_file, method, burnin, thin):
    labels = cluster_pyclone_trace(labels_file, method, burnin, thin)
    
    writer = csv.DictWriter(
                            open(cluster_file, 'w'),
                            ['mutation_id', 'cluster_id'],
                            delimiter='\t'
                            )
    
    writer.writeheader()
    
    for mutation_id, cluster_id in labels.items():
        out_row = {'mutation_id' : mutation_id, 'cluster_id' : int(cluster_id)}
        
        writer.writerow(out_row)  
        
def cluster_pyclone_trace(labels_file, method, burnin, thin):    
    trace = load_cluster_labels_trace(labels_file, burnin, thin)
    
    X = trace.values()
    
    if method == 'single_linkage':
        labels = cluster_with_hierachical(X, 'single')
    elif method == 'complete_linkage':
        labels = cluster_with_hierachical(X, 'complete')
    elif method == 'average_linkage':
        labels = cluster_with_hierachical(X, 'average')        
    elif method == 'affinity_propogation':
        labels = cluster_with_affinity_propogation(X)
    elif method == 'dbscan':
        labels = cluster_with_dbscan(X)
    elif method == 'spectral_clustering':
        labels = cluster_with_spectral_clustering(X)
    elif method == 'dynamic_tree_cut':
        labels = cluster_with_dynamic_tree_cut(X)
    elif method == 'mpear':
        labels = cluster_with_mpear(X)
    else:
        raise Exception("Clustering method {0} not recognised.".format(method))
    
    mutation_ids = trace.keys()
    
    return dict(zip(mutation_ids, labels))

def cluster_with_affinity_propogation(X):
    try:
        import sklearn.cluster as cluster
    except:
        raise Exception('''Clustering with affinity propogation requires the scikit-learn package. See http://scikit-learn.org/stable/.''')
    
    X = np.array(X)
    
    dist_mat = pdist(X, 'hamming')
    
    dist_mat = squareform(dist_mat)
    
    sim_mat = 1 - dist_mat
    
    ap = cluster.AffinityPropagation()
    
    ap.fit(sim_mat)
    
    return ap.labels_

def cluster_with_dbscan(X):
    try:
        import sklearn.cluster as cluster
    except:
        raise Exception('''Clustering with DBSCAN requires the scikit-learn package. See http://scikit-learn.org/stable/.''')
        
    X = np.array(X)
    
    dist_mat = pdist(X, 'hamming')
    
    dist_mat = squareform(dist_mat)
    
    dbscan = cluster.DBSCAN(metric='precomputed')
    
    dbscan.fit(dist_mat)
    
    labels = dbscan.labels_
    
    max_label = max(labels) + 1
    
    fixed_labels = []
    
    for label in labels:
        if label == -1:
            fixed_labels.append(max_label)
            
            max_label += 1
        else:
            fixed_labels.append(label)
    
    return fixed_labels

def cluster_with_dynamic_tree_cut(X):
    try:         
        import rpy2
    except:
        raise Exception('''Clustering with dynamic tree cut requires rpy2 package. See http://rpy.sourceforge.net/rpy2.html.''')
    
    import rpy2.robjects as robjects
        
    from rpy2.robjects.packages import importr    
    from rpy2.robjects.numpy2ri import numpy2ri
    
    robjects.conversion.py2ri = numpy2ri

    X = np.array(X)
    
    dist_mat = pdist(X, 'hamming')
    
    dist_mat = squareform(dist_mat)
    
    sim_mat = 1 - dist_mat
    
    importr('dynamicTreeCut')

    r = robjects.r
    
    robjects.globalenv['sim_mat'] = r.matrix(sim_mat, nrow=sim_mat.shape[0])
    
    r(
      '''
          dist.m <- function(d) dist(d, method="euclidean");
          hclust.m <- function(h) hclust(h, method="single");
          clusters <- cutreeDynamic(hclust.m(dist.m(sim_mat)), distM=sim_mat, minClusterSize=1, method="hybrid");
       '''
    )
    
    return np.array(r.clusters)

def cluster_with_hierachical(X, method):
    clusters = fclusterdata(X, 0.99, method=method, metric='hamming')
    
    return clusters

def cluster_with_spectral_clustering(X):
    try:
        import sklearn.cluster as cluster
    except:
        raise Exception('''Clustering with spectral clustering requires the scikit-learn package. See http://scikit-learn.org/stable/.''')
            
    X = np.array(X)
    
    dist_mat = pdist(X, 'hamming')
    
    dist_mat = squareform(dist_mat)
    
    sim_mat = 1 - dist_mat
    
    sc = cluster.SpectralClustering()
    
    sc.fit(sim_mat)
    
    return sc.labels_
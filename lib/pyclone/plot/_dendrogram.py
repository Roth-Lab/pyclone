'''
Created on 2012-08-05

@author: Andrew Roth
'''
from __future__ import division

from collections import defaultdict

import matplotlib.pyplot as plot
import numpy as np
import scipy.cluster.hierarchy as sch

class Dendrogram(object):
    def __init__(self, ax, method='single', cluster_threshold=None):
        self._ax = ax
        
        self.method = method
        
        self.cluster_threshold = cluster_threshold
        
        [i.set_linewidth(0) for i in self._ax.spines.itervalues()]
    
    def plot(self, genes, sim_mat):
        old_ax = plot.gca()
        
        plot.sca(self._ax)
        
        dist_mat = self._get_distance_matrix(sim_mat)
        
        linkage = sch.linkage(dist_mat, method=self.method)
        
        dendrogram = sch.dendrogram(linkage,
                                    orientation='right',
                                    labels=genes,
                                    count_sort=True,
                                    color_threshold=self.cluster_threshold)    
        
        self._ax.set_xticks([])
        
        plot.sca(old_ax)
        
        self.indices = dendrogram['leaves'][::-1]
    
    def get_clusters(self, genes, sim_mat):
        sim_mat = np.array(sim_mat)
        
        dist_mat = self._get_distance_matrix(sim_mat)
        
        linkage = sch.linkage(dist_mat, method=self.method)
        
        if self.cluster_threshold is None:
            cluster_threshold = 0.7 * np.max(linkage[:, 2])
        else:
            cluster_threshold = self.cluster_threshold  
        
        labels = sch.fcluster(linkage, cluster_threshold, 'distance')
    
        clusters = defaultdict(list)
            
        for i, gene in enumerate(genes):
            gene_label = labels[i]
        
            clusters[gene_label].append(gene)

        return clusters.values()
    
    def _get_distance_matrix(self, sim_mat):
        N = sim_mat.shape[0]
        
        dist_mat = 1 - sim_mat / np.amax(sim_mat, axis=1).reshape((1, N))
        
        return dist_mat
'''
Created on 2012-08-05

@author: Andrew Roth
'''
from __future__ import division

import matplotlib.pyplot as plot
import scipy.cluster.hierarchy as sch

class Dendrogram(object):
    def __init__(self, ax, method='average', cluster_threshold=0.99):
        self._ax = ax
        
        self.method = method
        
        self.cluster_threshold = cluster_threshold
        
        [i.set_linewidth(0) for i in self._ax.spines.itervalues()]
    
    def plot(self, genes, sim_mat):
        old_ax = plot.gca()
        
        plot.sca(self._ax)
        
        dist_mat = 1 - sim_mat
        
        linkage = sch.linkage(dist_mat, method=self.method)
        
        dendrogram = sch.dendrogram(linkage,
                                    orientation='right',
                                    labels=genes,
                                    count_sort=True,
                                    color_threshold=self.cluster_threshold)    
        
        self._ax.set_xticks([])
        
        plot.sca(old_ax)
        
        self.indices = dendrogram['leaves'][::-1]
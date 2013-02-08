'''
Created on 2012-08-05

@author: Andrew Roth
'''
from __future__ import division

from scipy.spatial.distance import pdist, squareform 

from ._colorbar import ColorBar
from ._dendrogram import Dendrogram
from ._heat_map import HeatMap

import matplotlib.pyplot as plot
import numpy as np

class SimilarityMatrixPlot(object):
    def __init__(self, labels_trace):
        self.trace = labels_trace
        
        self.genes = labels_trace.keys()

        self._init_plot_area()
    
    def plot(self):
        sim_mat = self._get_similarity_matrix()
        genes = np.array(self.genes)
        
        self.dendrogram.plot(genes, sim_mat)
        
        sim_mat = sim_mat[:, self.dendrogram.indices]
        sim_mat = sim_mat[self.dendrogram.indices, :][::-1]
        
        genes = genes[self.dendrogram.indices]
        
        self.heat_map.plot(genes, sim_mat)
        
        self.color_bar.plot(self.heat_map.palette)
        
        self.heat_map.x_labels = genes
        
        width = 0.5 * len(genes) + 3        
        height = 0.25 * len(genes) + 3
        
        self._fig.set_size_inches(width, height) 
        
    def save(self, file_name):
        self._fig.savefig(file_name, bbox_inches='tight')

    @property
    def title(self):
        return self.heat_map.title
    
    @title.setter
    def title(self, title):        
        self.heat_map.title = title    

    def _init_plot_area(self):
        self._fig = plot.figure()
        
        x_min = 0.05
        y_min = 0.15
        
        dendrogram_width = 0.1
        heat_map_width = 0.64
        color_bar_width = 0.02
        
        dend_to_hm = 0.06
        
        height = 0.7

        self._color_bar_ax = self._fig.add_axes([x_min + dendrogram_width + heat_map_width + dend_to_hm + 0.01,
                                                 y_min,
                                                 color_bar_width,
                                                 height])

        self._dendrogram_ax = self._fig.add_axes([x_min,
                                                  y_min,
                                                  dendrogram_width,
                                                  height])
        
        self._heat_map_ax = self._fig.add_axes([x_min + dendrogram_width + dend_to_hm,
                                                y_min,
                                                heat_map_width,
                                                height])
    
        self.color_bar = ColorBar(self._color_bar_ax)
        
        self.dendrogram = Dendrogram(self._dendrogram_ax)
        
        self.heat_map = HeatMap(self._heat_map_ax)            

    def _get_similarity_matrix(self):
        X = self.trace.values()
        
        dist_mat = pdist(X, 'hamming')
        
        dist_mat = squareform(dist_mat)
        
        sim_mat = 1 - dist_mat
        
        return sim_mat

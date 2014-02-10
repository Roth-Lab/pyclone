'''
Created on 2012-08-05

@author: Andrew Roth
'''
from __future__ import division

from collections import OrderedDict

from .densities import PosteriorDensity

import matplotlib.cm as cm
import matplotlib.pyplot as plot
import numpy as np

class CellularFrequencyPlot(object):
    def __init__(self, frequencies_trace, clusters=None, cmap=cm.Set1):
        self.trace = frequencies_trace        
        
        self.clusters = clusters        
        
        self.cmap = cmap

        self._init_pdfs()
                
        self._init_plot_area()
        
        self._init_clusters()
                
    def plot(self, sort_clusters=True, sort_genes=True):
        i = 0
        
        genes = []
        
        if sort_clusters:
            clusters = self._sort_clusters()
        else:
            clusters = self.clusters
        
        if sort_genes:
            for  cluster, color in zip(clusters, self.colors):
                for gene in cluster:
                    pos = self.positions[i]
                    
                    self.pdfs[gene].plot(self._ax, pos, self.width, color=color)
                    
                    i += 1
                    
                    genes.append(gene)
        else:
            for gene in self.pdfs:
                pos = self.positions[i]
                                                
                color = self.gene_color[gene]
                
                self.pdfs[gene].plot(self._ax, pos, self.width, color=color)
                
                i += 1
            
                genes.append(gene)
        
        self._fix_axes(genes)
        
        self._ax.set_xlabel('Cellular Frequency')
        
        self._ax.set_xlim(-0.01, 1.01)

    def save(self, file_name):
        self._fig.tight_layout()
        
        self._fig.savefig(file_name, bbox_inches='tight')
        
    @property
    def title(self):
        return self._title
    
    @title.setter
    def title(self, title):
        self._title = title
        
        self._ax.set_title(title)       

    def _init_pdfs(self):
        self.pdfs = OrderedDict()

        for gene in self.trace:
            self.pdfs[gene] = PosteriorDensity(self.trace[gene])
    
    def _init_plot_area(self):
        self._fig = plot.figure()
        
        self._ax = self._fig.add_subplot(1, 1, 1)
        
        self.num_loci = len(self.trace)        
        
        self.positions = range(1, self.num_loci + 1)
        
        self.distance = max(self.positions) - min(self.positions)
        
        self.width = min(0.15 * max(self.distance, 1.0), 0.5)        
        
        width = 11
          
        height = 0.25 * self.num_loci + 3
        
        self._fig.set_size_inches(width, height)
    
    def _init_clusters(self):
        self.gene_color = {}
        
        if self.clusters == None:
            self.clusters = [self.trace.keys(), ]
                    
        self.num_clusters = len(self.clusters)
        
        self.colors = [self.cmap(i / 10, 1) for i in range(self.num_clusters)]
        
        for i, cluster in enumerate(self.clusters):
            for gene in cluster:
                self.gene_color[gene] = self.colors[i]        
    
    def _fix_axes(self, genes):
        new_limits = (min(self.positions) - 0.5, max(self.positions) + 0.5)
        
        self._ax.set_ylim(new_limits)
        
        self._ax.set_yticks(self.positions)
        
        self._ax.set_yticklabels(genes, fontsize=8)

        box = self._ax.get_position()
    
        self._ax.set_position([box.x0 + box.width * 0.2, box.y0, box.width * 0.8, box.height])
                
    def _sort_clusters(self):
        sort_objective = np.mean
        
        cluster_values = [[] for _ in self.clusters]
        
        for i, cluster in enumerate(self.clusters):
            for gene in cluster:
                cluster_values[i].extend(self.trace[gene])
        
        # Sort clusters in order
        sort_values = []
        
        for value in cluster_values:
            sort_values.append(sort_objective(value))
        
        new_clusters = []
        
        for _, cluster in sorted(zip(sort_values, self.clusters)):
            new_clusters.append(cluster)
        
        # Sort within clusters
        for i, cluster in enumerate(new_clusters):
            sort_values = [sort_objective(self.trace[gene]) for gene in cluster]
            
            new_clusters[i] = []
            
            for _, gene in sorted(zip(sort_values, cluster)):
                new_clusters[i].append(gene)
            
        return new_clusters   

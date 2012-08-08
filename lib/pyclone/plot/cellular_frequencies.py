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
    def __init__(self, values, clusters, cmap=cm.Set1, burnin=1000, thin=10):
        self._fig = plot.figure()
        
        self._ax = self._fig.add_subplot(1, 1, 1)
        
        self.pdfs = OrderedDict()     
        self.values = OrderedDict()
        
        for gene in values:
            self.values[gene] = values[gene][burnin::thin]
            self.pdfs[gene] = PosteriorDensity(self.values[gene])
                
        self.clusters = clusters
        
        self.num_clusters = len(clusters)
        
        self.num_cols = len(values)
        
        self.positions = range(1, self.num_cols + 1)
        
        self.distance = max(self.positions) - min(self.positions)
        
        self.width = min(0.15 * max(self.distance, 1.0), 0.5)
        
        self.colors = [cmap(i / 10, 1) for i in range(self.num_clusters)]
        
        width = 11        
        height = 0.25 * len(values) + 3
        
        self._fig.set_size_inches(width, height)
        
        self.gene_color = {}
        
        for i, cluster in enumerate(self.clusters):
            for gene in cluster:
                self.gene_color[gene] = self.colors[i] 
        
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
        
        new_limits = (min(self.positions) - 0.5, max(self.positions) + 0.5)
        
        self._ax.set_ylim(new_limits)
        
        self._ax.set_yticks(self.positions)

#        short_gene_names = [x.split('_')[0] for x in genes]
        
        self._ax.set_yticklabels(genes, fontsize=8)

        box = self._ax.get_position()
    
        self._ax.set_position([box.x0 + box.width * 0.2, box.y0, box.width * 0.8, box.height])
        
        self._ax.set_xlabel('Cellular Frequency')

    def save(self, file_name):
        self._fig.savefig(file_name, bbox_inches='tight')
        
    @property
    def title(self):
        return self._title
    
    @title.setter
    def title(self, title):
        self._title = title
        
        self._ax.set_title(title)       
        
    def _sort_clusters(self):
        sort_objective = np.mean
        
        cluster_values = [[] for _ in self.clusters]
        
        for i, cluster in enumerate(self.clusters):
            for gene in cluster:
                cluster_values[i].extend(self.values[gene])
        
        # Sort clusters in order
        sort_values = []
        
        for value in cluster_values:
            sort_values.append(sort_objective(value))
        
        new_clusters = []
        
        for _, cluster in sorted(zip(sort_values, self.clusters)):
            new_clusters.append(cluster)
        
        # Sort within clusters
        for i, cluster in enumerate(new_clusters):
            sort_values = [sort_objective(self.values[gene]) for gene in cluster]
            
            new_clusters[i] = []
            
            for _, gene in sorted(zip(sort_values, cluster)):
                new_clusters[i].append(gene)
            
        return new_clusters   

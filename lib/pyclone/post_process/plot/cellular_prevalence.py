'''
Created on 2012-08-05

@author: Andrew Roth
'''
from __future__ import division

from collections import OrderedDict

import matplotlib.pyplot as pp

from pyclone.post_process import cluster_pyclone_trace
from pyclone.post_process.utils import load_cellular_frequencies_trace

from .densities import PosteriorDensity
from .utils import get_clusters_color_map, setup_axes, setup_plot

def plot_cellular_prevalence_posteriors(labels_trace_file, prevalence_trace_file, plot_file, burnin, thin):  
    labels = cluster_pyclone_trace(labels_trace_file, burnin, thin)
    
    labels = labels.set_index('mutation_id')
    
    labels = labels.sort_values(by='cluster_id')
    
    labels = labels['cluster_id']
  
    prevalence_trace = load_cellular_frequencies_trace(prevalence_trace_file, burnin, thin)
  
    plotter = CellularFrequencyPlot(labels, prevalence_trace)

    plotter.plot()
    
    plotter.save(plot_file)

class CellularFrequencyPlot(object):
    def __init__(self, labels, prevalence_trace):
        self.clusters = labels.unique()
        
        self.labels = labels
        
        self.loci = labels.index
        
        self.trace = prevalence_trace        

        self.num_loci = self.trace.shape[1]
        
        self._init_locus_colors()

        self._init_pdfs()
                
        self._init_plot_area()
                
    def plot(self):
        for i, locus in enumerate(self.loci):
            pos = self.positions[i]
            
            self.pdfs[locus].plot(self._ax, pos, self.width, color=self.locus_colors[locus])
   
        self._fix_axes()
        
        self._ax.set_xlabel('Cellular prevalence')
        
        self._ax.set_ylabel('Loci')
        
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

        for locus in self.trace:
            self.pdfs[locus] = PosteriorDensity(self.trace[locus])
    
    def _init_plot_area(self):
        setup_plot()
        
        self._fig = pp.figure()
        
        self._ax = self._fig.add_subplot(1, 1, 1)
        
        setup_axes(self._ax)
      
        self.positions = range(1, self.num_loci + 1)
        
        self.distance = max(self.positions) - min(self.positions)
        
        self.width = min(0.15 * max(self.distance, 1.0), 0.5)        
        
        width = 8
          
        height = 0.25 * self.num_loci + 3
        
        self._fig.set_size_inches(width, height)
    
    def _init_locus_colors(self):            
        color_map = get_clusters_color_map(self.labels)
        
        self.locus_colors = self.labels.map(color_map)

    def _fix_axes(self):
        new_limits = (min(self.positions) - 0.5, max(self.positions) + 0.5)
        
        self._ax.set_ylim(new_limits)
        
        self._ax.set_yticks(self.positions)
        
        self._ax.set_yticklabels(self.loci, fontsize=8)

        box = self._ax.get_position()
    
        self._ax.set_position([box.x0 + box.width * 0.2, box.y0, box.width * 0.8, box.height])

'''
Created on 2012-08-05

@author: Andrew Roth
'''
from __future__ import division

from .densities import PosteriorDensity

import brewer2mpl
import matplotlib.pyplot as plot
import pandas as pd

# Setup color paletter
bmap = brewer2mpl.get_map('Set1', 'qualitative', 9)

palette = bmap.mpl_colors

class CellularFrequencyPlot(object):
    def __init__(self, frequencies_trace):
        self.trace = frequencies_trace        

        self._init_plot_area()
                
    def plot(self):
        mean_cellular_prevalence = self.trace.mean(axis=1)
        
        mean_cellular_prevalence = pd.DataFrame({'mean_cellular_prevalence' : mean_cellular_prevalence}, 
                                                index=self.trace.index)
        
        mean_cellular_prevalence.sort(columns='mean_cellular_prevalence', inplace=True)

        row_ids = mean_cellular_prevalence.index
        
        for i, row_id in enumerate(row_ids):
            pos = self.positions[i]
            
            pdf = PosteriorDensity(self.trace.loc[row_id])
                                            
            pdf.plot(self._ax, pos, self.width, color=palette[0])
            
        self._fix_axes(row_ids)
        
        self._ax.set_xlabel('Cellular Prevalence')
        
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
    
    def _fix_axes(self, genes):
        new_limits = (min(self.positions) - 0.5, max(self.positions) + 0.5)
        
        self._ax.set_ylim(new_limits)
        
        self._ax.set_yticks(self.positions)
        
        self._ax.set_yticklabels(genes, fontsize=8)

        box = self._ax.get_position()
    
        self._ax.set_position([box.x0 + box.width * 0.2, box.y0, box.width * 0.8, box.height])
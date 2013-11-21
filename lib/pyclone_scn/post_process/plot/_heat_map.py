'''
Created on 2012-08-05

@author: Andrew Roth
'''
import matplotlib.cm as cm

class HeatMap(object):
    def __init__(self, ax):
        self._ax = ax
        
        self.palette = None
        
        self._ax.tick_params(length=0)
        
    def plot(self, genes, sim_mat, color_map='Blues'):
        self.palette = self._ax.pcolor(sim_mat, cmap=cm.get_cmap(color_map), edgecolors='k', linewidth='1.5')
        
        self.x_labels = genes
        self.y_labels = []
        
        self._fix_ticks(genes)  
        
    @property
    def title(self):
        return self._title
    
    @title.setter
    def title(self, title):
        self._title = title
        
        self._ax.set_title(title)
        
    @property
    def x_labels(self):
        return self._x_labels
    
    @x_labels.setter
    def x_labels(self, labels):
        self._x_labels = labels
        
        self._ax.set_xticklabels(labels, size=10, rotation=90)

    @property
    def y_labels(self):
        return self._y_labels
    
    @y_labels.setter
    def y_labels(self, labels):
        self._y_labels = labels
        
        self._ax.set_yticklabels(labels, size=8)
    
    def _fix_ticks(self, labels):
        positions = [x + .5 for x in range(0, len(labels))]
    
        newlimits = min(positions) - 0.5, max(positions) + 0.5

        self._ax.set_xticks(positions)        
        self._ax.set_xlim(newlimits)
    
        self._ax.set_yticks(positions)        
        self._ax.set_ylim(newlimits)        
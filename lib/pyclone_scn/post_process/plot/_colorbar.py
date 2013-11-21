'''
Created on 2012-08-05

@author: Andrew Roth
'''
import matplotlib.pyplot as plot

class ColorBar(object):
    def __init__(self, ax):
        self._ax = ax
        
        self._label = None
        
        self._color_bar = None
    
    def plot(self, palette):
        self._color_bar = plot.colorbar(palette, cax=self._ax)
        
    @property
    def label(self):
        return self._label
    
    @label.setter
    def label(self, label):
        self._label = label
        
        ticks = [x[0] for x in label]
        tick_labels = [x[1] for x in label]
        
        self._color_bar.set_ticks(ticks)
        self._color_bar.set_ticklabels(tick_labels)  
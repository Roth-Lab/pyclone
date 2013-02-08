'''
Created on 2012-08-05

@author: Andrew Roth
'''
from __future__ import division

import numpy as np
import scipy.stats as stats

class PosteriorDensity(object):
    def __init__(self, sample_values, support=(0, 1), grid_size=1000):
        '''
        Fit a gaussian kernel density estimator to a set sample points for plotting.
        
        Args:
            sample_values : (list) Values from MCMC sampler.
        
        Kwargs:
            support : (tuple) The support of the density.
            grid_size : (int) Number of points to evaluate density at.
        '''
        sample_values = np.array(sample_values)
        
        sample_values += np.random.uniform(0, 1e-6, len(sample_values))
        
        self.kde = stats.gaussian_kde(sample_values)
        
        self.x = np.linspace(support[0], support[1], grid_size)
        
        self.y = self.kde.evaluate(self.x)
        
    @property
    def pdf(self):
        return self.x, self.y
    
    def plot(self, ax, position, width, color='y', alpha=0.8):
        x = self.x
        
        y = width * self.y / self.y.max()
        
        ax.fill_between(x, position, position + y, facecolor=color, alpha=alpha)

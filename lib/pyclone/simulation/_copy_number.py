'''
Created on 2012-06-14

@author: Andrew Roth
'''
from __future__ import division

from numpy.random import dirichlet

from pyclone.utils import discrete_rvs

class CopyNumberSimulator(object):
    '''
    Object for simulating copy numbers.
    '''
    def __init__(self, min_cn, max_cn, mix_weight_priors=None):
        '''
        Args:
            min_cn (int): Minimum copy number to simulate.
            max_cn (int): Maximum copy number to simulate.
        
        Kwargs:
            mix_weight_priors (list of ints): List of Dirichlet prior pseudo counts for the copy number classes. Entries
                                              start at min_cn and go to max_cn.
        '''
        self.min_cn = min_cn
        self.max_cn = max_cn
        
        self.num_classes = max_cn - min_cn + 1
        
        if mix_weight_priors is None:
            mix_weight_priors = [2., ] * self.num_classes
        else:
            if len(mix_weight_priors) != self.num_classes:
                raise Exception("Length of priors list should match the number of copy number classes.")
        
        self.mix_weights = dirichlet(mix_weight_priors)
        
    def draw_cn(self):
        '''
        Draw a copy number.
        '''
        cn = discrete_rvs(self.mix_weights)
        
        # Shift cn from 0 based to min_cn
        cn = cn + self.min_cn
        
        return cn
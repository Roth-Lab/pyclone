'''
Created on 2011-12-29

@author: Andrew Roth
'''
from __future__ import division

from math import log
from random import random

from clonal_estimation.model import CollapsedLikelihood

class CollapsedSampler(object):
    def run(self, data_point, max_iters=100000, thin=1, burnin=0):
        results = {'phi' : []}

        likelihood = CollapsedLikelihood(data_point)
            
        phi = random()
        
        for i in range(max_iters):
            phi = self._update_phi(likelihood, phi)
            
            if i % thin == 0 and i >= burnin:
                results['phi'].append(phi)
                print i, phi
        
        return results

    def _update_phi(self, likelihood, old_phi):   
        new_phi = random()
        
        numerator = likelihood.compute_log_likelihood(new_phi)
        
        denominator = likelihood.compute_log_likelihood(old_phi)
        
        log_ratio = numerator - denominator
        
        u = random()
        
        if log_ratio >= log(u):
            return new_phi
        else:
            return old_phi
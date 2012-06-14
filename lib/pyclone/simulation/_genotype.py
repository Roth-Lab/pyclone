'''
Created on 2012-06-14

@author: Andrew Roth
'''
from __future__ import division

from numpy.random import dirichlet

from pyclone.utils import discrete_rvs

class GenotypeSimulator(object):
    '''
    Object for simulating a genotype like AB or AAB.
    '''
    def __init__(self, cn_simulator, mix_weight_priors=None):
        '''
        Args:
            cn_simulator (:class:CopyNumberSimulator): Copy number simulator used to generate copy number for genotypes.
        
        Kwargs:
            mix_weight_priors (list of lists): Each entry should be a list of Dirichlet prior pseudo counts for the
            genotypes assosciated with a copy number. Entries should go from min_cn to max_cn specified in cn_simulator
            object.
        '''
        self._cn_simulator = cn_simulator
        
        self.min_cn = cn_simulator.min_cn
        self.max_cn = cn_simulator.max_cn
        
        if mix_weight_priors is None:
            mix_weight_priors = {}
            
            for cn in range(self.min_cn, self.max_cn + 1):
                mix_weight_priors[cn] = [2., ] * cn        
        else:
            if len(mix_weight_priors) != self._cn_simulator.num_classes:
                raise Exception('''Length of mix-weight priors list should be the same as the number of copy number
                                classes.''')
        
        self.mix_weights = {}
        
        for cn in mix_weight_priors:
            self.mix_weights[cn] = dirichlet(mix_weight_priors[cn])
    
    def draw_genotype(self):
        '''
        Draw a genotype such as AB or AAB.
        '''
        cn = self._cn_simulator.draw_cn()
                
        num_ref_alleles = discrete_rvs(self.mix_weights[cn])
        
        genotype = "A" * num_ref_alleles + "B" * (cn - num_ref_alleles)
        
        return genotype 
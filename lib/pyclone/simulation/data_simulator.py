'''
Created on 2012-01-29

@author: Andrew Roth
'''
from __future__ import division

from collections import defaultdict

from pyclone.utils import discrete_rvs, binomial_rvs, poisson_rvs

import random

class ClonalDataSimulator(object):
    def __init__(self,
                 num_clones,
                 min_cn,
                 max_cn,
                 mean_depth,
                 random_cn_mix_weights=False,
                 random_genotype_mix_weights=False):
        
        self._cn_simulator = CopyNumberSimulator(min_cn, max_cn, random_cn_mix_weights)
        self._genotype_simulator = GenotypeSimulator(min_cn, max_cn, random_genotype_mix_weights)
         
        self.mean_depth = mean_depth
        
        self.mu_ref = 0.999
        
        self.clone_mix_weights = [1 / num_clones] * num_clones
        self.clone_frequencies = [random.random() for _ in range(num_clones)]
    
    def draw_data_point(self):
        clone_label = discrete_rvs(self.clone_mix_weights)
        
        clone_freq = self.clone_frequencies[clone_label]
        
        cn = self._cn_simulator.draw_cn()
        
        genotype = self._genotype_simulator.draw_genotype(cn)
        
        mu_var = self._genotype_simulator.get_reference_frequency(cn, genotype)
        
        mu = (1 - clone_freq) * self.mu_ref + clone_freq * mu_var
        
        d = poisson_rvs(self.mean_depth)
        
        a = binomial_rvs(mu, d)
        
        return a, d, cn, genotype, clone_freq

class CopyNumberSimulator(object):
    def __init__(self, min_cn, max_cn, random_mix_weights=False):
        self.min_cn = min_cn
        self.max_cn = max_cn
        
        self.num_classes = max_cn - min_cn + 1
        
        # Use either random mixture weights or uniform.
        if random_mix_weights:
            self.mix_weights = [random.random() for _ in range(self.num_classes)]
            self.mix_weights = [x / sum(self.mix_weights) for x in self.mix_weights]
        else:
            self.mix_weights = [1 / self.num_classes] * self.num_classes
        
    def draw_cn(self):
        cn = discrete_rvs(self.mix_weights)
        
        # Shift cn from 0 based to min_cn
        cn = cn + self.min_cn
        
        return cn

class GenotypeSimulator(object):
    def __init__(self, min_cn, max_cn, random_mix_weights=False):
        self.min_cn = min_cn
        self.max_cn = max_cn
        
        self.num_cn_classes = max_cn - min_cn + 1
        
        self._init_mix_weights(random_mix_weights)
        
        self._init_reference_allele_frequencies()
    
    def _init_mix_weights(self, random_mix_weights):
        self.mix_weights = defaultdict(list)
        
        for cn in range(self.min_cn, self.max_cn + 1):
            num_genotypes = cn + 1
            
            if random_mix_weights:
                self.mix_weights[cn] = [random.random() for _ in range(num_genotypes)]                 
            else:
                self.mix_weights[cn] = [1 / num_genotypes] * num_genotypes
            
            # Assign AA class 0 mix-weight
            self.mix_weights[cn][-1] = 0
            
            # Normalise mix-weights
            self.mix_weights[cn] = [ x / sum(self.mix_weights[cn]) for x in self.mix_weights[cn]]
    
    def _init_reference_allele_frequencies(self):
        self.reference_allele_frequencies = defaultdict(list)
        
        for cn in range(self.min_cn, self.max_cn + 1):
            num_genotypes = cn + 1
            
            self.reference_allele_frequencies[cn] = [x / cn for x in range(num_genotypes)]
    
    def draw_genotype(self, cn):
        genotype = discrete_rvs(self.mix_weights[cn])
        
        return genotype
    
    def get_reference_frequency(self, cn, genotype):
        return self.reference_allele_frequencies[cn][genotype]
                             
if __name__ == "__main__":
    simulator = ClonalDataSimulator(10, 2, 2, 1000)
    
    for i in range(100):
        print simulator.draw_data_point()

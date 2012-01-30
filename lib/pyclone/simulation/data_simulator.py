'''
Created on 2012-01-29

@author: Andrew Roth
'''
from __future__ import division

from collections import defaultdict

from pyclone.utils import discrete_rvs

from scipy.stats import binom as binomial, poisson

import random

class SimulatorFactory(object):
    @staticmethod
    def get_simple_simulator(min_cn, max_cn, num_clones, mean_depth):
        cn_sim = CopyNumberSimulator(min_cn, max_cn)
        
        genotype_sim = GenotypeSimulator(cn_sim)
        
        clone_sim = CloneSimulator(num_clones)
        
        sim = ClonalDataSimulator(clone_sim, genotype_sim, mean_depth)
        
        return sim
        
class ClonalDataSimulator(object):
    def __init__(self, clone_simulator, genotype_simulator, mean_depth, mu_ref=0.999):        
        self._clone_simulator = clone_simulator
        self._genotype_simulator = genotype_simulator
         
        self.mean_depth = mean_depth
        
        self.mu_ref = mu_ref

    def draw_data_point(self):       
        d = poisson.rvs(self.mean_depth)
        
        clone_freq = self._clone_simulator.draw_clonal_frequency()
                
        d_var = binomial.rvs(d, clone_freq)
        d_ref = d - d_var
        
        mu_ref = self.mu_ref
        a_ref = binomial.rvs(d_ref, mu_ref)
        
        mu_var, genotype = self._genotype_simulator.draw_mu()
        a_var = binomial.rvs(d_var, mu_var)
        
        a = a_var + a_ref
        
        return a, d, genotype, clone_freq
    
    @property
    def clone_frequencies(self):
        return self._clone_simulator.frequencies
    
    @clone_frequencies.setter
    def clone_frequencies(self, value):
        self._clone_simulator.frequencies = value

class CloneSimulator(object):
    def __init__(self, num_clones, random_mix_weights=False):
        self.num_clones = num_clones
        
        if random_mix_weights:
            self.mix_weights = [random.random() for _ in range(self.num_clones)]
            self.mix_weights = [x / sum(self.mix_weights) for x in self.mix_weights]
        else:
            self.mix_weights = [1 / self.num_clones] * self.num_clones
        
        self.frequencies = [random.random() for _ in range(num_clones)]
    
    def draw_clonal_frequency(self):
        population_id = discrete_rvs(self.mix_weights)
        
        return self.frequencies[population_id]

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
    def __init__(self, cn_simulator, random_mix_weights=False, error_rate=0.001):
        self._cn_simulator = cn_simulator
        
        self.min_cn = cn_simulator.min_cn
        self.max_cn = cn_simulator.max_cn
        
        self.error_rate = error_rate
        
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
        self.mu = defaultdict(list)
        
        for cn in range(self.min_cn, self.max_cn + 1):
            num_genotypes = cn + 1
            
            self.mu[cn] = [x / cn for x in range(num_genotypes)]
            self.mu[cn][0] = self.error_rate
            self.mu[cn][-1] = 1 - self.error_rate
    
    def draw_mu(self):
        cn = self._cn_simulator.draw_cn()    
        
        genotype = discrete_rvs(self.mix_weights[cn])
        
        return self.mu[cn][genotype], (cn, genotype)
                             
if __name__ == "__main__":
    sim = SimulatorFactory.get_simple_simulator(2, 2, 1, 1000)
    
    for i in range(100):
        print sim.draw_data_point()

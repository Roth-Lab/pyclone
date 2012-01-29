'''
Created on 2012-01-29

@author: Andrew Roth
'''
from __future__ import division

from pyclone.utils import discrete_rvs, binomial_rvs, poisson_rvs

import random

class ClonalDataSimulator(object):
    def __init__(self, params):
        self.params = params
    
    def draw_data_point(self):
        clone_label = discrete_rvs(self.params.clone_mix_weights)
        
        clone_freq = self.params.clone_frequencies[clone_label]
        
        cn = discrete_rvs(self.params.cn_mix_weights)
        
        genotype_mix_weights = self.params.genotype_mix_weights[cn]
        
        genotype = discrete_rvs(genotype_mix_weights)
        
        mu = (1 - clone_freq) * self.params.mu_ref + clone_freq * self.params.mu_var[cn][genotype]
        
        d = poisson_rvs(self.params.mean_depth)
        
        a = binomial_rvs(mu, d)
        
        return a, d, cn, genotype

class ClonalDataSimulatorParamters(object):
    def __init__(self, num_clones, max_cn, mean_depth):
        self.mean_depth = mean_depth
        
        self.mu_ref = 0.999
        
        self.clone_mix_weights = [1 / num_clones] * num_clones
        self.clone_frequencies = [random.random() for _ in range(num_clones)]
        
        self.cn_mix_weights = [0, ] + [1 / max_cn] * max_cn
        
        self.genotype_mix_weights = []
        self.mu_var = []
        
        for cn in range(max_cn + 1):
            if cn == 0:
                self.genotype_mix_weights.append([1, ])
                self.mu_var.append([0, ])
            else:
                num_genotypes = cn + 1
                mix_weights = []
                mu_var = []
                for g in range(num_genotypes):
                    mix_weights.append(1 / (num_genotypes))
                    mu_var.append(g / (num_genotypes))
                
                self.genotype_mix_weights.append(mix_weights)
                self.mu_var.append(mu_var)
        
                    
if __name__ == "__main__":
    params = ClonalDataSimulatorParamters(10, 8, 1000)
    simulator = ClonalDataSimulator(params)
    
    print simulator.draw_data_point()

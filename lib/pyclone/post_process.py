'''
Created on 2012-02-08

@author: Andrew Roth
'''
from __future__ import division

from collections import defaultdict

class PostProcessor(object):
    def __init__(self, data):       
        self.genes = data['genes']
        
        self._db = data
        
        self._results = data['trace']
        
    @property
    def alpha(self):
        return self._results['alpha']
    
    @property
    def cellular_frequencies(self):
        '''
        Returns a dictionary with keys genes, and values posterior samples of cellular frequencies.
        '''
        phi = defaultdict(list)
        
        labels = self.labels
        
        for gene in labels:
            for label, sample in zip(labels[gene], self._results['phi']):
                phi[gene].append(sample[label])
        
        return phi

    @property
    def labels(self):
        '''
        Returns a dict with keys genes, and values the class label of the genes for each MCMC sample.
        '''
        labels = defaultdict(list)
        
        for sample in self._results['labels']:
            for gene, label in zip(self.genes, sample):
                labels[gene].append(label)
        
        return labels

    @property
    def num_components(self):
        '''
        Returns a list of the number of components used in by each MCMC sample.
        '''
        labels = self._results['labels']
        
        num_components = []
        
        for sample in labels:
            num_components.append(len(set(sample)))
        
        return num_components
    
    @property
    def num_iterations(self):
        '''
        Returns the number of MCMC iterations.
        '''
        return self._db['sampler'].num_iters
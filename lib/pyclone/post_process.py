'''
Created on 2012-02-08

@author: Andrew Roth
'''
from __future__ import division

from collections import defaultdict, Counter
from math import ceil, floor, log

from pyclone.utils import histogram, log_sum_exp
from pyclone.model import BinomialLikelihood

class DpSamplerPostProcessor(object):
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
    def likelihood(self):
        '''
        Return a dictionary of likelihood objects, one object for each gene.
        '''
        likelihoods = {}
        
        for gene, data_point in zip(self._db['genes'], self._db['data']):
            likelihoods[gene] = BinomialLikelihood(data_point)
        
        return likelihoods

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
    def num_samples(self):
        return self._db['sampler'].num_iters

    def get_similarity_matrix(self, burnin=0, thin=1):
        '''
        Gets the posterior similarity matrix. The i,j entry is the number of times gene i and j where in the same
        cluster.
        '''            
        n = len(self.genes)
        
        sim_mat = [[0] * n for _ in range(n)] 
        
        labels = self.labels
        
        for i in range(n):
            for j in range(i, n):
                gene_1 = self.genes[i]
                gene_2 = self.genes[j]
                
                sim_mat[i][j] = self._get_gene_similarity(labels[gene_1][burnin::thin], labels[gene_2][burnin::thin])
                sim_mat[j][i] = sim_mat[i][j]
        
        return sim_mat

    def get_alpha_posteriors(self, burnin=0, thin=1):        
        alpha = self.alpha[burnin::thin]
        
        min_value = floor(min(alpha))
        max_value = ceil(max(alpha))
        
        num_bins = int(max_value - min_value)
        
        return histogram(alpha, num_bins, min_value=min_value, max_value=max_value, normalise=True)

    def get_cellular_frequency_posteriors(self, num_bins=100, burnin=0, thin=1):
        posteriors = {}
        
        cellular_frequencies = self.cellular_frequencies
        
        for gene in cellular_frequencies:
            posteriors[gene] = histogram(cellular_frequencies[gene][burnin::thin],
                                         num_bins,
                                         min_value=0,
                                         max_value=1,
                                         normalise=True)
        
        return posteriors
    
    def get_num_component_posteriors(self, burnin=0, thin=1):
        num_components = self.num_components[burnin::thin]
        
        min_value = floor(min(num_components)) - 0.5       
        max_value = ceil(max(num_components)) + 0.5
        
        num_bins = int(max_value - min_value)
                
        return histogram(num_components, num_bins, min_value=min_value, max_value=max_value, normalise=True)        
    
    def get_similarity_posteriors(self, burnin=0, thin=1):
        sim_mat = self.get_similarity_matrix(burnin=burnin, thin=thin)
        
        for i, row in enumerate(sim_mat):
            sim_mat[i] = [x / max(row) for x in row]
        
        return sim_mat

    def _get_gene_similarity(self, labels_1, labels_2):
        '''
        Computes how frequently items in two lists match.
        '''
        similarity = 0
        
        for l1, l2 in zip(labels_1, labels_2):
            if l1 == l2:
                similarity += 1
        
        return similarity
    
    @property
    def map_clusters(self):
        likelihoods = self.likelihood
                
        frequencies = self.cellular_frequencies
        
        labels = self.labels
        
        partitions = self._results['labels']
        
        alpha = self.alpha
        
        max_ll = float('-inf')
        best_sample = None
        
        for i in range(self.num_samples):
            ll = []
            
            partition_size = Counter(partitions[i])
            
            for gene in likelihoods:
                phi = frequencies[gene][i]
                gene_ll = likelihoods[gene].compute_log_likelihood(phi)
                
                gene_partition = labels[gene][i]
                
                weight = log(partition_size[gene_partition] + alpha[i] - 1)
                
                ll.append(gene_ll + weight)
            
            ll = log_sum_exp(ll)
            
            if ll > max_ll:
                max_ll = ll
                best_sample = i
        
        clusters = defaultdict(list)
        
        for gene in frequencies:
            clusters[frequencies[gene][best_sample]].append(gene)
        
        return clusters
            

if __name__ == "__main__":
    post_processor = DpSamplerPostProcessor("../../examples/test.pickle")
    
    print post_processor.similarity_matrix
    print post_processor.num_components


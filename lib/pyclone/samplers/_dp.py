'''
Created on 2011-12-15

@author: Andrew Roth
'''
from __future__ import division

from math import exp, log
from random import betavariate as beta_rvs, gammavariate as gamma_rvs, random, shuffle

from pyclone.utils import bernoulli_rvs, discrete_rvs, log_space_normalise

class DirichletProcessSampler(object):
    def __init__(self, data, m=2, concentration=None):
               
        self._clusters = Clusters(data)
                                                
        self._seat_sampler = LabelUpdater(m)            
        
        self._dish_sampler = FrequencyUpdater()
        
        if concentration is None:
            self._concentration_sampler = ConcentrationUpdater(1e-3, 1e-3)
        
            self._update_concentration = True
        else:
            self._update_concentration = False
            
            self._seat_sampler.concentration_parameter = concentration
        
    def sample(self, num_iters=1000, burnin=0, thin=1):
        print self._clusters.num_members
        
        results = {'alpha' : [], 'labels' : [], 'phi' : []}
        
        for i in range(num_iters):
            self._update_phi()
            self._update_labels()
            
            if self._update_concentration:
                self._update_concentration_parameters()
            
            if i % thin == 0 and i >= burnin:
                print i, self._clusters, self._seat_sampler.concentration_parameter
                
                results['alpha'].append(self._seat_sampler.concentration_parameter)
                results['labels'].append(self._clusters.labels)
                results['phi'].append(self._clusters.values)
        
        return results
    
    def  _update_phi(self):
        self._dish_sampler.update_clusters(self._clusters)
        
    def _update_labels(self):
        self._seat_sampler.update_clusters(self._clusters)

    def _update_concentration_parameters(self):
        conc_param = self._concentration_sampler.draw(self._seat_sampler.concentration_parameter,
                                                      self._clusters.num_clusters,
                                                      self._clusters.num_members)
        
        self._seat_sampler.concentration_parameter = conc_param

class LabelUpdater(object):
    '''
    Chinese restaurant process sampling using algorithm 8 from Neal.
    '''
    def __init__(self, m, concentration_parameter=1):
        self.m = m
        
        self.concentration_parameter = concentration_parameter

    def update_clusters(self, clusters):
        # Shuffle the order of position to reseat.
        members = range(clusters.num_members)
        shuffle(members)
        
        for i in members:
            self._reseat_member(i, clusters)
    
    def _reseat_member(self, member, clusters):
        member_value = clusters.pop_member(member)
        
        counts = self._get_counts(clusters)
        values = self._get_values(clusters, member_value)

        data = clusters.data[member]
        
        label = self._get_label(data, counts, values)
        
        if label < clusters.num_clusters:
            clusters.push_member(member, label)
        else:                        
            clusters.add_cluster(member, values[label])
    
    def _get_counts(self, clusters):
        counts = clusters.counts
        
        counts.extend([self.concentration_parameter / self.m] * self.m)
        
        return counts

    def _get_values(self, clusters, member_value):
        values = clusters.values[:]
                
        values.extend(self._get_temp_phi(clusters, member_value))
        
        return values

    def _get_temp_phi(self, clusters, member_value):
        temp_phi = []
        
        if clusters.is_cluster_empty(member_value):
            temp_phi.append(member_value)
        
        while len(temp_phi) < self.m:
            temp_phi.append(random())
        
        return temp_phi
    
    def _get_label(self, data, counts, values):
        p = self._get_crp_probs(data, counts, values)

        label = discrete_rvs(p)
        
        return label
    
    def _get_crp_probs(self, data, counts, values):
        h = len(counts)
        
        log_p = []
        
        for k in range(h):
            log_p.append(log(counts[k]) + data.compute_log_likelihood(values[k]))

        log_p = log_space_normalise(log_p)
        
        p = [exp(x) for x in log_p]
        
        return p
    
class FrequencyUpdater(object):
    '''
    Class for sampling the clonal frequencies (dishes) using metropolis step.
    '''    
    def update_clusters(self, clusters):
        old_values = clusters.values
        
        new_values = []
    
        for k, old_phi in enumerate(old_values):            
            new_phi = random()
            
            cluster_data_points = clusters.get_cluster_data_points(k)
            
            if self._accept_new_phi(new_phi, old_phi, cluster_data_points):
                new_values.append(new_phi)
            else:
                new_values.append(old_phi)
        
        clusters.values = new_values
    
    def _accept_new_phi(self, new_phi, old_phi, cluster_data):
        '''
        Do MH step.
        '''        
        numerator = 0.
        
        denominator = 0.
    
        for data_point in cluster_data:
            numerator += data_point.compute_log_likelihood(new_phi)
       
            denominator += data_point.compute_log_likelihood(old_phi)
            
        log_ratio = numerator - denominator

        u = random()
        
        if log_ratio >= log(u):
            return True
        else:
            return False

class ConcentrationUpdater(object):
    '''
    Class for sampling alpha parameters in Dirichlet process, using a gamma prior with paramters a,b.
    '''
    def __init__(self, a, b):
        self.a = a
        self.b = b
            
    def draw(self, old_alpha, k, n):
        eta = self._get_eta(old_alpha, n)
        
        pi = self._get_pi(eta, k, n)
        
        # Add a small positive number to avoid errors
        alpha = self._get_alpha(eta, pi, k) + 1e-10
        
        return alpha
        
    def _get_eta(self, old_alpha, n):
        return beta_rvs(old_alpha + 1, n)
    
    def _get_pi(self, eta, k, n):
        a = self.a
        b = self.b
    
        x = (a + k - 1) / (n * (b - log(eta)))
        
        return x / (1 + x)
    
    def _get_alpha(self, eta, pi, k):
        a = self.a
        b = self.b
        
        label = bernoulli_rvs(pi)
        
        scale = 1 / (b - log(eta))
                
        if label == 0:
            return gamma_rvs(a + k, scale)
        elif label == 1:
            return gamma_rvs(a + k - 1, scale)
        else:
            raise Exception("Invalid label drawn.")

#=======================================================================================================================
# 
#=======================================================================================================================
class Clusters(object):
    def __init__(self, data):
        self.data = data
        
        self.num_members = len(data)
        self.num_clusters = 0

        self.tables = []
        self.values = []
        
        self._init_clusters()
    
    def __str__(self):
        return str(self.num_clusters)
    
    @property
    def labels(self):
        labels = [-1] * self.num_members
        
        for k, table in enumerate(self.tables):
            for member in table:
                labels[member] = k
        
        return labels

    @property
    def counts(self):
        counts = [len(table) for table in self.tables]
        
        return counts
        
    def add_cluster(self, member, value):
        self.values.append(value)
        self.tables.append([member, ])
    
        self.num_clusters = len(self.tables)
    
    def pop_member(self, member):
        value = -1
        
        for k, table in enumerate(self.tables):
            if member in table:
                table.remove(member)
                value = self.values[k]
            
                if len(table) == 0:
                    self._remove_cluster(k)
                
                break
        
        if value == -1:
            raise Exception('Clusters.pop_member : Member {0} not found.'.format(member))
        
        return value
    
    def push_member(self, member, cluster_id):
        self.tables[cluster_id].append(member)
     
    def is_cluster_empty(self, value):
        if value in self.values:
            return 0
        else:
            return 1
        
    def get_cluster_data_points(self, cluster_id):
        return [self.data[i] for i in self.tables[cluster_id]]
    
    
    def _init_clusters(self):        
        for i in range(self.num_members):
            value = random()
            
            self.add_cluster(i, value)

    def _remove_cluster(self, index):        
        self.tables.pop(index)
        self.values.pop(index)
        
        self.num_clusters = len(self.tables)

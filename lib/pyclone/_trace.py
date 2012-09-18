'''
Created on 2012-05-10

@author: Andrew
'''
from collections import defaultdict

import os

import pyclone.zshelve as shelve

class TraceDB(object):
    def __init__(self, file_prefix, mode='r', max_cache_size=1000):
        if mode in ['a', 'r'] and not os.path.exists(file_prefix):
            raise Exception("{0} does not exists.".format(file_prefix))
        elif mode == 'w':
            if os.path.exists(file_prefix):
                raise Exception("{0} exists, cannot overwrite.".format(file_prefix))
            if not os.path.exists(os.path.dirname(os.path.abspath(file_prefix))):
                raise Exception("Folder {0} does not exist to create pyclone file in.".format(os.path.dirname(file_prefix)))
        
        self.mode = mode
        
        self._load_db(file_prefix)

        self._cache_size = 0
            
        self._max_cache_size = max_cache_size
    
    def _load_db(self, file_prefix):
        '''
        Load the shelve db object if it exists, otherwise initialise.
        '''
        if self.mode == 'r':
            self._db = shelve.open(file_prefix, writeback=False)
        else:
            self._db = shelve.open(file_prefix, writeback=True)
        
        # Check if file exists, if not initialise
        if 'trace' not in self._db:                
            self._db['trace'] = {'alpha' : [], 'labels' : [], 'phi' : []}
            
            self._db.sync()
                
    def __getitem__(self, key):
        return self._db[key]
    
    def __setitem__(self, key, value):
        if self.mode == 'r':
            raise Exception('AnalysisDB cannot be edited in read only mode.')
             
        self._db[key] = value 
        
    def update_trace(self, state):
        for parameter in self._db['trace']:
            self._db['trace'][parameter].append(state[parameter])
        
        self._cache_size += 1
        
        if self._cache_size >= self._max_cache_size:
            self._db.sync()
            self._cache_size = 0
    
    def close(self):
        self._db.close() 
    
    def sync(self):
        self._db.sync()

class TracePostProcessor(object):
    def __init__(self, file_name, burnin, thin):
        trace_db = TraceDB(file_name, mode='r')
        
        self.genes = trace_db['genes'][:]
            
        self.alpha = trace_db['trace']['alpha'][burnin::thin]
        
        self.dishes = trace_db['trace']['phi'][burnin::thin]
        
        self.partitions = trace_db['trace']['labels'][burnin::thin]
        
        self.num_iterations = trace_db['sampler'].num_iters
        
        trace_db.close()

    @property
    def cellular_frequencies(self):
        '''
        Returns a dictionary with keys genes, and values posterior samples of cellular frequencies.
        '''
        phi = defaultdict(list)
        
        labels = self.labels
        
        dishes = self.dishes
        
        for gene in labels:
            for label, dish_sample in zip(labels[gene], dishes):
                phi[gene].append(dish_sample[label])
        
        return phi

    @property
    def labels(self):
        '''
        Returns a dict with keys genes, and values the class label of the genes for each MCMC sample.
        '''
        labels = defaultdict(list)

        for partition in self.partitions:
            for gene, cluster_id in zip(self.genes, partition):
                labels[gene].append(cluster_id)
        
        return labels

    @property
    def num_components(self):
        '''
        Returns a list of the number of components used in by each MCMC sample.
        '''        
        num_components = []
        
        for partition in self.partitions:
            num_components.append(len(set(partition)))
        
        return num_components
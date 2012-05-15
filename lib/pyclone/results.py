'''
Created on 2012-05-10

@author: Andrew
'''
import shelve

class SamplerResults(object):
    def __init__(self, file_name, mode='r', max_cache_size=1000):
        self._db = shelve.open(file_name, writeback=True)
        
        if mode == 'w':                
            self._db['trace'] = {'alpha' : [], 'labels' : [], 'phi' : []}
            
            self._db.sync()
            
        self._cache_size = 0
            
        self._max_cache_size = max_cache_size
    
    def __getitem__(self, key):
        return self._db[key]
    
    def __setitem__(self, key, value):        
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

'''
Created on 2012-05-10

@author: Andrew
'''
import os
import platform
import shelve

class SamplerResults(object):
    def __init__(self, file_prefix, max_cache_size=1000):
        self._load_db(file_prefix)

        self._cache_size = 0
            
        self._max_cache_size = max_cache_size
    
    def _load_db(self, file_prefix):
        '''
        Load the shelve db object if it exists, otherwise initialise.
        ''' 
        # Workaround for PyPY shelve implementation adding .db and Cpython not
        if platform.python_implementation() == "PyPy":
            file_name = file_prefix 
        else:
            file_name += ".db"
        
        self._db = shelve.open(file_name, writeback=True)
        
        # Check if file exists, if not initialise
        if not os.path.exists(file_name):                
            self._db['trace'] = {'alpha' : [], 'labels' : [], 'phi' : []}
            
            self._db.sync()
        else:
            # If file exists and trace is not in db raise error.
            if 'trace' not in self._db:
                raise Exception('{0} exists but does not appear to be a valid PyClone SamplerResults object.'.format(file_name))
                
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

'''
Created on 2013-04-23

@author: Andrew Roth
'''
import bz2
import csv
import os

from pyclone.utils import make_directory

class DiskTrace(object):
    def __init__(self, trace_dir, sample_ids, mutation_ids, attribute_map, precision=False):
        self.trace_dir = trace_dir
        
        self.sample_ids = sample_ids
        
        self.mutation_ids = mutation_ids
        
        self.attribute_map = attribute_map
        
        self.update_precision = precision 
    
    def close(self):
        self.alpha_writer.close()
        
        self.labels_writer.close()
            
        for writer in self.cellular_frequency_writers.values():
            writer.close()
            
        if self.update_precision:
            self.precision_writer.close()
    
    def open(self):
        make_directory(self.trace_dir)
        
        self.alpha_writer = ConcentrationParameterWriter(self.trace_dir)
        
        self.labels_writer = LabelsWriter(self.trace_dir, self.mutation_ids)
        
        self.cellular_frequency_writers = {}
        
        for sample_id in self.sample_ids:
            self.cellular_frequency_writers[sample_id] = CellularFrequenciesWriter(self.trace_dir, 
                                                                                   sample_id, 
                                                                                   self.mutation_ids)
        
        if self.update_precision:
            self.precision_writer = PrecisionWriter(self.trace_dir)
    
    def update(self, state):
        self.alpha_writer.write_row([state['alpha'], ])
        
        self.labels_writer.write_row(state['labels'])
        
        for sample_id in self.sample_ids:
            row = []
            
            for param in state['params']:
                attr = self.attribute_map['cellular_frequencies']
                
                row.append(getattr(param[sample_id], attr))
            
            self.cellular_frequency_writers[sample_id].write_row(row)
            
        if self.update_precision:
            self.precision_writer.write_row([state['global_params'].x])        

class ConcentrationParameterWriter(object):
    def __init__(self, trace_dir):
        self.file_name = os.path.join(trace_dir, 'alpha.tsv.bz2')
    
        self.file_handle = bz2.BZ2File(self.file_name, 'w')
        
        self.writer = csv.writer(self.file_handle, delimiter='\t')
        
        self.param_id = 'alpha'
    
    def close(self):
        self.file_handle.close()
    
    def write_row(self, row):
        self.writer.writerow(row)

class CellularFrequenciesWriter(object):
    def __init__(self, trace_dir, sample_id, mutation_ids):
        self.file_name = os.path.join(trace_dir, '{0}.cellular_prevalence.tsv.bz2'.format(sample_id))
    
        self.file_handle = bz2.BZ2File(self.file_name, 'w')
        
        self.writer = csv.writer(self.file_handle, delimiter='\t')
        
        self.writer.writerow(mutation_ids)
        
        self.param_id = (sample_id, 'cellular_frequencies')
    
    def close(self):
        self.file_handle.close()
    
    def write_row(self, row):
        self.writer.writerow(row)

class LabelsWriter(object):
    def __init__(self, trace_dir, mutation_ids):
        self.file_name = os.path.join(trace_dir, 'labels.tsv.bz2')
    
        self.file_handle = bz2.BZ2File(self.file_name, 'w')
        
        self.writer = csv.writer(self.file_handle, delimiter='\t')
        
        self.writer.writerow(mutation_ids)
        
        self.param_id = 'labels'
    
    def close(self):
        self.file_handle.close()
    
    def write_row(self, row):
        self.writer.writerow(row)

class PrecisionWriter(object):
    def __init__(self, trace_dir):
        self.file_name = os.path.join(trace_dir, 'precision.tsv.bz2')
    
        self.file_handle = bz2.BZ2File(self.file_name, 'w')
        
        self.writer = csv.writer(self.file_handle, delimiter='\t')
        
        self.param_id = 'precision'
    
    def close(self):
        self.file_handle.close()
    
    def write_row(self, row):
        self.writer.writerow(row)        
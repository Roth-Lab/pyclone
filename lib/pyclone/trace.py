'''
Created on 2013-04-23

@author: Andrew Roth
'''
import bz2
import csv

from pyclone.utils import make_directory

import pyclone.paths as paths

class DiskTrace(object):
    def __init__(self, config_file, mutation_ids, attribute_map, precision=False):
        self.config_file = config_file
    
        self.sample_ids = paths.get_sample_ids(config_file)
        
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
        make_directory(paths.get_trace_dir(self.config_file))
        
        self.alpha_writer = ConcentrationParameterWriter(paths.get_concentration_trace_file(self.config_file))
        
        self.labels_writer = LabelsWriter(
                                paths.get_labels_trace_file(self.config_file), 
                                self.mutation_ids
                            )
        
        self.cellular_frequency_writers = {}
        
        for sample_id, file_name in paths.get_cellular_prevalence_trace_files(self.config_file).items():
            self.cellular_frequency_writers[sample_id] = CellularFrequenciesWriter(
                                                            file_name,
                                                            sample_id,
                                                            self.mutation_ids
                                                        )
        
        if self.update_precision:
            self.precision_writer = PrecisionWriter(paths.get_precision_trace_file(self.config_file))
    
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
    def __init__(self, file_name):
        self.file_name = file_name
    
        self.file_handle = bz2.BZ2File(self.file_name, 'w')
        
        self.writer = csv.writer(self.file_handle, delimiter='\t')
        
        self.param_id = 'alpha'
    
    def close(self):
        self.file_handle.close()
    
    def write_row(self, row):
        self.writer.writerow(row)

class CellularFrequenciesWriter(object):
    def __init__(self, file_name, mutation_ids, sample_id):
        self.file_name = file_name
    
        self.file_handle = bz2.BZ2File(self.file_name, 'w')
        
        self.writer = csv.writer(self.file_handle, delimiter='\t')
        
        self.writer.writerow(mutation_ids)
        
        self.param_id = (sample_id, 'cellular_frequencies')
    
    def close(self):
        self.file_handle.close()
    
    def write_row(self, row):
        self.writer.writerow(row)

class LabelsWriter(object):
    def __init__(self, file_name, mutation_ids):
        self.file_name = file_name
    
        self.file_handle = bz2.BZ2File(self.file_name, 'w')
        
        self.writer = csv.writer(self.file_handle, delimiter='\t')
        
        self.writer.writerow(mutation_ids)
        
        self.param_id = 'labels'
    
    def close(self):
        self.file_handle.close()
    
    def write_row(self, row):
        self.writer.writerow(row)

class PrecisionWriter(object):
    def __init__(self, file_name):
        self.file_name = file_name
    
        self.file_handle = bz2.BZ2File(self.file_name, 'w')
        
        self.writer = csv.writer(self.file_handle, delimiter='\t')
        
        self.param_id = 'precision'
    
    def close(self):
        self.file_handle.close()
    
    def write_row(self, row):
        self.writer.writerow(row)        

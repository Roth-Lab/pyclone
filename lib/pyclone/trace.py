'''
Created on 2013-04-23

@author: Andrew Roth
'''
import bz2
import csv
import os
import pandas as pd

from pyclone.utils import make_directory

class HDF5Trace(object):
    def __init__(self, file_name, mutation_ids, sample_ids, store_precision=False):
        self._file_name = file_name
        
        self.mutation_ids = mutation_ids
        
        self.sample_ids = sample_ids
        
        self.store_precision = store_precision
    
    def __enter__(self):
        self.open()
        
        return self
    
    def __exit__(self, type, value, traceback):
        self.close()
    
    def open(self):
        self.trace = pd.HDFStore(self._file_name, 'w')#, complevel=9, complib='blosc')
        
        self._init_labels_table()
        
        for sample_id in self.sample_ids:
            self._init_cellular_prevalence_table(sample_id)
        
        self._init_parameter_table('alpha')
        
        if self.store_precision:
            self._init_parameters_table('precision')
    
    def close(self):
        self.trace.close()
    
    def update(self, state):
        for sample_id in self.sample_ids:
            values = []
             
            for param in state['params']:
                values.append(param[sample_id].x)
             
            self._update_cellular_prevalence_table(sample_id, values)
         
        self._update_labels_table(state['labels'])
        
        self._update_parameter_table('alpha', state['alpha'])
         
        if self.store_precision:
            self._update_parameter_table('precision', state['global_params'].x)
    
    def _init_labels_table(self):
        df = pd.DataFrame(columns=self.mutation_ids, dtype='int32')
        
        self.trace.append('labels', df)
    
    def _init_cellular_prevalence_table(self, sample_id):
        df = pd.DataFrame(columns=self.mutation_ids, dtype='float64')
        
        table_name = 'cellular_prevalences/{0}'.format(sample_id)
        
        self.trace.append(table_name, df)
    
    def _init_parameter_table(self, parameter_name):
        df = pd.DataFrame(columns=['value'])
        
        df['value'] = df['value'].astype('float64')
        
        table_name = 'parameters/{0}'.format(parameter_name)
        
        self.trace.append(table_name, df)
    
    def _update_cellular_prevalence_table(self, sample_id, values):
        table_name = 'cellular_prevalences/{0}'.format(sample_id)
        
        df = pd.DataFrame(values, index=self.mutation_ids)
        
        df = df.T
        
        self.trace.append(table_name, df)
    
    def _update_labels_table(self, values):
        table_name = 'labels'
        
        df = pd.DataFrame(values, index=self.mutation_ids)
        
        df = df.T
        
        self.trace.append(table_name, df)
    
    def _update_parameter_table(self, parameter_name, value):
        table_name = 'parameters/{0}'.format(parameter_name)
        
        df = pd.DataFrame([value], index=['value'])
        
        df = df.T
        
        self.trace.append(table_name, df)

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
        self.file_name = os.path.join(trace_dir, '{0}.cellular_frequencies.tsv.bz2'.format(sample_id))
    
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
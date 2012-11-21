'''
Created on 2012-05-10

@author: Andrew
'''
import bz2
import csv
import os

class TraceDB(object):
    def __init__(self, out_dir, mutations):
        if os.path.exists(out_dir):
            raise Exception("{0} exists, cannot overwrite.".format(out_dir))
    
        if not os.path.exists(os.path.dirname(os.path.abspath(out_dir))):
            raise Exception("Folder {0} does not exist to create pyclone file in.".format(os.path.dirname(out_dir)))
        
        os.makedirs(out_dir)
  
        self._open_files(out_dir)
        
        self._frequencies_writer.writerow(mutations)
        
        self._labels_writer.writerow(mutations)

    def _open_files(self, out_dir):
        '''
        Load the shelve db object if it exists, otherwise initialise.
        '''
        mode = 'w'
        
        self._alpha_file = bz2.BZ2File(os.path.join(out_dir, 'alpha.tsv.bz2'), mode)
        
        self._frequencies_file = bz2.BZ2File(os.path.join(out_dir, 'cellular_frequencies.tsv.bz2'), mode)
        
        self._labels_file = bz2.BZ2File(os.path.join(out_dir, 'labels.tsv.bz2'), mode)
        
        self._phi_file = bz2.BZ2File(os.path.join(out_dir, 'phi.tsv.bz2'), mode)        

        self._alpha_writer = csv.writer(self._alpha_file, delimiter='\t')
            
        self._frequencies_writer = csv.writer(self._frequencies_file, delimiter='\t')
        
        self._labels_writer = csv.writer(self._labels_file, delimiter='\t')
            
        self._phi_writer = csv.writer(self._phi_file, delimiter='\t')
        
    def update_trace(self, state):
        self._alpha_writer.writerow([state['alpha'], ])
        
        self._frequencies_writer.writerow(state['cellular_frequencies'])
        
        self._labels_writer.writerow(state['labels'])
        
        self._phi_writer.writerow(state['phi'])

    def close(self):
        self._alpha_file.close()
        
        self._frequencies_file.close()
        
        self._labels_file.close()
        
        self._phi_file.close()

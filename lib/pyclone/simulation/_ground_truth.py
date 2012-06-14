'''
Created on 2012-06-14

@author: Andrew Roth
'''
from __future__ import division

import csv

from collections import OrderedDict

from pyclone.simulation import CopyNumberSimulator, GenotypeSimulator, MixtureModel, ChineseRestaurantProcess

class ClonalMutationDataSet(object):
    '''
    Object representing a data-set comprised of mutations at sub-clonal cellular frequencies.
    
    Currently only supports uniformly random genotype and cellular frequency distributions.
    '''    
    def generate_from_cluster_model(self, num_mutations, num_clones, min_cn, max_cn):
        '''
        Args:
            num_mutations (int): Number of mutations in data set.
            num_clones (int): Number of clonal populations.
            min_cn (int): Minimum copy number of genotype.
            max_cn (int): Maximum copy number of genotype.
        '''
        cn_simulator = CopyNumberSimulator(min_cn, max_cn)
    
        genotype_simulator = GenotypeSimulator(cn_simulator)
    
        cellular_freq_simulator = MixtureModel(num_clones)
        
        self.data = OrderedDict()
        
        for i in range(num_mutations):
            gene = "gene_{0}".format(i)
            
            genotype = genotype_simulator.draw_genotype()
            
            cellular_freq = cellular_freq_simulator.draw_cellular_frequency()
            
            self.data[gene] = ClonalMuationDataPoint(genotype, cellular_freq)

    def generate_from_crp_model(self, num_mutations, alpha, min_cn, max_cn):
        '''
        Args:
            num_mutations (int): Number of mutations in data set.
            alpha (float): Concentration parameter in CRP, lower values will result in fewer clusters.
            min_cn (int): Minimum copy number of genotype.
            max_cn (int): Maximum copy number of genotype.
        '''
        cn_simulator = CopyNumberSimulator(min_cn, max_cn)
    
        genotype_simulator = GenotypeSimulator(cn_simulator)
    
        cellular_freq_simulator = ChineseRestaurantProcess(alpha, num_mutations)
        
        self.data = OrderedDict()
        
        for i in range(num_mutations):
            gene = "gene_{0}".format(i)
            
            genotype = genotype_simulator.draw_genotype()
            
            cellular_freq = cellular_freq_simulator.draw_cellular_frequency()
            
            self.data[gene] = ClonalMuationDataPoint(genotype, cellular_freq)            
    
    def load_from_file(self, file_name):
        reader = csv.DictReader(open(file_name), delimiter='\t')
        
        self.data = OrderedDict()
        
        for row in reader:
            gene = row['gene']
            
            genotype = row['genotype']
            
            cellular_freq = float(row['cellular_frequency'])
            
            data_point = ClonalMuationDataPoint(genotype, cellular_freq)
            
            self.data[gene] = data_point
                
    def save_to_file(self, file_name):
        '''
        Save data-set to file. Format is tab delimited with header. Fields are gene, genotype, cellular_frequency.
        
        Args:
            file_name (str): Path to where the file will be save.
        '''
        fields = ['gene', 'genotype', 'cellular_frequency']
        
        writer = csv.DictWriter(open(file_name, 'w'), fields, delimiter='\t')
        writer.writeheader()
        
        for gene, data_point in self.data.items():
            out_row = {}
            
            out_row['gene'] = gene
            
            out_row['genotype'] = data_point.genotype
            
            out_row['cellular_frequency'] = data_point.cellular_freq
            
            writer.writerow(out_row)

class ClonalMuationDataPoint(object):
    '''
    Object representing a mutation at a sub-clonal frequency.
    '''
    def __init__(self, genotype, cellular_freq):
        '''
        Args:
            genotype (str): String representation of genotype i.e. AAB
            cellular_freq (float): Cellular frequency of variant population at position.
        '''
        self.genotype = genotype        
        
        self.cellular_freq = cellular_freq
        
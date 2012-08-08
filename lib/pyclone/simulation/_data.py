'''
Created on 2012-06-14

@author: Andrew Roth
'''
from __future__ import division

import csv

from collections import defaultdict, OrderedDict
from numpy.random import binomial, poisson

class PycloneSyntheticDataSet(object):
    '''
    Object representing a sample from a pyclone data set. The data is count data, LOH, and copy number data with some
    noise like would be generated from NGS sequencing experiment.
    ''' 
    def __init__(self, clonal_population_data, mean_depth, error_rate=0.001):
        '''
        Args:
           clonal_population_data (:class: ClonalMutationDataSet): The true population which will be sampled from.
           mean_depth (int): Average depth of coverage. This is the mean parameter for Poisson depth sampling.
        
        Kwargs:
            error_rate (float): Error rate of sequencing technology i.e. how frequently an A is called as a B or vice
                                versa.
        '''        
        self.data = OrderedDict()
        
        for gene, data_point in clonal_population_data.data.items():
            # Draw depth above 0
            while True: 
                depth = poisson(mean_depth)
                
                if depth > 0:
                    break
            
            self.data[gene] = PycloneSyntheticDataPoint(data_point, depth, error_rate)
    
    def save_to_file(self, file_name, prior_model='vague'):
        '''
        Save data-set to file. Format Pyclone compatible tab delimited with header. Fields are gene, a, d, mu_r,
        delta_r, mu_v, delta_v.
        
        Args:
            file_name (str): Path to where the file will be save.

        Kwargs:
            prior_model (str): Model to be used to assign priors valid values ['vague', 'informative','diploid_no_loh',
                               'diploid_loh'].
        '''        
        fields = ['gene', 'a', 'd', 'mu_r', 'delta_r', 'mu_v', 'delta_v', 'genotype', 'cellular_frequency']
        
        writer = csv.DictWriter(open(file_name, 'w'), fields, delimiter='\t')
        
        writer.writeheader()
        
        for gene, data_point in self.data.items():            
            out_row = {}
            
            out_row['gene'] = gene
            
            out_row['a'] = data_point.a
            out_row['d'] = data_point.d
            
            out_row['mu_r'] = self._list_to_comma_separated_string(data_point.mu_ref)
            out_row['delta_r'] = self._list_to_comma_separated_string(data_point.delta_ref)
            
            mu_var, delta_var = data_point.get_variant_priors(prior_model) 
            
            out_row['mu_v'] = self._list_to_comma_separated_string(mu_var)
            out_row['delta_v'] = self._list_to_comma_separated_string(delta_var)
            
            out_row['genotype'] = data_point.genotype
            out_row['cellular_frequency'] = data_point.cellular_frequency
            
            writer.writerow(out_row)

    def _list_to_comma_separated_string(self, x):
        return ",".join([str(x_i) for x_i in x])

class PycloneSyntheticDataPoint(object):
    def __init__(self, clonal_muatation_data_point, depth, error_rate):
        '''        
        Args:
            depth (int): Depth of coverage for this position.
        
            clonal_muatation_data_point (:class: ClonalMuationDataPoint): Data point to sample counts, genotype
                                                                          estimates for.
                                                                          
            error_rate (float): Error rate of sequencing technology i.e. how frequently an A is called as a B or vice
                                versa.            
        '''
        self.d = depth
        
        self._data_point = clonal_muatation_data_point
        
        self._error_rate = error_rate
        
        self.min_expected_cn = 1
        self.max_expected_cn = 6
        
        self._set_counts()
        
        self._set_reference_priors()
    
    def get_variant_priors(self, prior_model):
        '''
        prior_model (str): Model to be used to assign priors valid values ['vague', 'informative','diploid_no_loh',
                           'diploid_loh'].
        '''
        if prior_model == 'vague':
            priors = self._get_vague_priors()
        elif prior_model == 'informative':
            priors =  self._get_informative_priors()
        elif prior_model == 'diploid_no_loh':
            priors =  self._get_diploid_no_loh_priors()
        elif prior_model == 'diploid_loh':
            priors = self._get_diploid_loh_priors()
        else:
            raise Exception("Prior model `{0}` not recognised.".format(prior_model))
        
        return priors
    
    def _set_counts(self):
        mu_var = self._get_ref_allele_sampling_frequency(self._data_point.genotype)
        
        a = self._draw_ref_counts(mu_var, self._data_point.cellular_freq)
        
        self.a = a

    def _get_ref_allele_sampling_frequency(self, genotype):
        '''
        Get the sampling frequency associated with a genotype.
        
        Args:
            genotype (str): Genotype to get sampling frequency for in form AAB.
        ''' 
        cn = len(genotype)
        
        num_ref_alleles = genotype.count('A')
        
        if num_ref_alleles == 0:
            sampling_freq = self._error_rate
        elif num_ref_alleles == cn:
            raise Exception('Genotype {0} not allowed in variant population.'.format(self.genotype))
        else:
            sampling_freq = num_ref_alleles / cn
        
        return sampling_freq

    def _draw_ref_counts(self, mu_var, cellular_freq):
        '''
        Draw the number of observed ref counts from this position.
        
        Args:
            mu_var (float): The sampling frequency associated with the genotype of the variant population.
            cellular_freq (float): Proportion of variant cells in the population.
        '''        
        d = self.d
        
        # Draw number of reads from ref and variant populations.
        d_var = binomial(d, cellular_freq)        
        d_ref = d - d_var
                
        mu_ref = 1 - self._error_rate
                
        if d_ref == 0:
            a_ref = 0
        else:
            a_ref = binomial(d_ref, mu_ref)
          
        if d_var == 0:
            a_var = 0
        else:
            a_var = binomial(d_var, mu_var)
        
        a = a_var + a_ref
        
        return a

    def _set_reference_priors(self):
        self.mu_ref = [1 - self._error_rate, ]
        self.delta_ref = [1, ]

    def _get_vague_priors(self):
        '''
        Convert copy number and LOH calls to a prior.
        
        Frequencies compatible with the predicted cn get the highest priors, with lower priors as the CN deviates from the
        predicted cn.
        
        If the site is LOH the homozygous reference state is assigned to be five times as likely. Otherwise the heterozygous
        states are five times as likely.
        '''
        pseudo_counts = defaultdict(int)
        
        true_cn = len(self._data_point.genotype)
        true_num_ref_alleles = self._data_point.genotype.count('A')

        if true_num_ref_alleles == 0:
            is_loh = True
        else:
            is_loh = False        

        # Use copy number information
        for cn in range(self.min_expected_cn, self.max_expected_cn + 1):
            for num_ref_alleles in range(cn):
                if num_ref_alleles == 0:
                    mu_v = self._error_rate
                else:
                    mu_v = num_ref_alleles / cn
                
                cn_error = abs(true_cn - cn)
                
                if cn_error == 0:
                    pseudo_counts[mu_v] += 10
                elif cn_error == 1:
                    pseudo_counts[mu_v] += 7.5
                elif cn_error == 2:
                    pseudo_counts[mu_v] += 5
                elif cn_error == 3:
                    pseudo_counts[mu_v] += 2.5
                else:
                    pseudo_counts[mu_v] += 1
        
        # Use LOH information
        for mu_v in pseudo_counts:
            if mu_v == self._error_rate and is_loh == 1:
                pseudo_counts[mu_v] *= 5
            elif mu_v != self._error_rate and is_loh == 0:
                pseudo_counts[mu_v] *= 5
        
        mu_var = sorted(pseudo_counts.keys())        
        delta_var = [pseudo_counts[x] for x in mu_var]
        
        return mu_var, delta_var
    
    def _get_informative_priors(self):
        true_cn = len(self._data_point.genotype)
        true_num_ref_alleles = self._data_point.genotype.count('A')
                
        if true_num_ref_alleles == 0:
            is_loh = True
        else:
            is_loh = False
        
        priors = defaultdict(int)
        
        for num_ref_alleles in range(true_cn):
            if num_ref_alleles == 0:
                mu_var = self._error_rate
            else:
                mu_var = num_ref_alleles / true_cn
            
            priors[mu_var] = 1
        
        if is_loh:
            priors[self._error_rate] = 10                
        elif not is_loh:
            priors[self._error_rate] = 1 / 10
        
        mu_var = sorted(priors.keys())
        delta_var = [priors[x] for x in mu_var]
        
        return mu_var, delta_var
        
    
    def _get_diploid_no_loh_priors(self):
        mu_var = [self._error_rate, 0.5]
        delta_var = [1, 1]
        
        return mu_var, delta_var
    
    def _get_diploid_loh_priors(self):
        mu_var = [self._error_rate, 0.5]
        
        true_num_ref_alleles = self._data_point.genotype.count('A')
                
        if true_num_ref_alleles == 0:
            is_loh = True
        else:
            is_loh = False
        
        if is_loh:
            delta_var = [10, 1]
        else:
            delta_var = [1 / 10, 1]
        
        return mu_var, delta_var
    
    @property
    def genotype(self):
        return self._data_point.genotype
    
    @property
    def cellular_frequency(self):
        return self._data_point.cellular_freq

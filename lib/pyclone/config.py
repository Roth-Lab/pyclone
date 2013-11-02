'''
Created on 2013-02-12

@author: Andrew Roth
'''
from __future__ import division

import csv

def load_mutations_from_csv(file_name):
    mutations = []
    
    reader = csv.DictReader(open(file_name), delimiter='\t')

    for row in reader:
        mutation_id = row['mutation_id']
        
        ref_counts = int(row['ref_counts'])
        
        var_counts = int(row['var_counts'])
        
        normal_cn = int(row['normal_cn'])
        
        major_cn = int(row['major_cn'])
        
        minor_cn = int(row['minor_cn'])
        
        cn_prevalence = float(row['cn_prevalence'])

        mutation = get_mutation(mutation_id,
                                ref_counts,
                                var_counts,
                                normal_cn,
                                major_cn,
                                minor_cn,
                                cn_prevalence)

        mutations.append(mutation)    
    
    return mutations

def get_mutation(mutation_id, ref_counts, var_counts, normal_cn, major_cn, minor_cn, cn_prevalence):    
    mutation = Mutation(mutation_id, ref_counts, var_counts, normal_cn, major_cn, minor_cn, cn_prevalence)

    return mutation

#=======================================================================================================================
# Helper classes
#=======================================================================================================================
class Mutation(object):
    def __init__(self, mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn, cn_prevalence):
        self.id = mutation_id
        
        self.ref_counts = ref_counts
        
        self.var_counts = var_counts
        
        self.normal_cn = normal_cn
        
        self.major_cn = major_cn
        
        self.minor_cn = minor_cn
        
        self.total_cn = major_cn + minor_cn
        
        self.cn_prevalence = cn_prevalence

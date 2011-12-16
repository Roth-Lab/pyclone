from __future__ import division

from clonal_estimation.independent_model import DataPoint, CollapsedLikelihood
from clonal_estimation.dp_model import DirichletProcessSampler
from math import exp

import csv
import multiprocessing
import os

def main(in_file_name, plot_folder):
    reader = csv.DictReader(open(in_file_name), delimiter='\t')
    
#    data = [run_gene(reader.next()) for _ in range(5)]
    data = [run_gene(row) for row in reader]

#    p = multiprocessing.Pool()
#    data = p.map(run_gene, [x for x in reader])    
#    p.close()
    
    sampler = DirichletProcessSampler(data, concentration=1)
    
    sampler.sample()
    
def run_gene(row, max_depth=100):
    gene = row['gene']
    
    print gene
    
    a = int(row['ref_counts'])
    b = int(row['var_counts'])
    
    d = a + b
    
    if d > max_depth:
        a = round((a / d) * max_depth)
        d = max_depth
        
        print a, d
        
    
    cn = int(row['cn'])
    loh = int(row['loh'])
    
    mu, pi_r, pi_v = get_genotype_priors(cn, loh)
    
    data_point = DataPoint(a, d, pi_r, pi_v, mu)
    
    n = 1000
    likelihood = CollapsedLikelihood(data_point)
    
    return likelihood

def get_genotype_priors(cn, loh):
    eps = 1e-3
    
    mu = []
    
    for i in range(cn + 1):
        if i == 0:
            mu.append(eps)
        elif i == cn:
            mu.append(1 - eps)
        else:
            mu.append(i / cn)
    
    pi_r = [0.] * len(mu)
    pi_r[-1] = 1
    
    pi_v = [0.] * len(mu)
    
    for i in range(len(mu) - 1):
        if loh:
            if i == 0:
                pi_v[i] = 10
            else:
                pi_v[i] = 1
        else:
            if i == 0:
                pi_v[i] = 1
            else:
                pi_v[i] = 10
    
    pi_v = [x / sum(pi_v) for x in pi_v]
    
    return mu, pi_r, pi_v        
        

if __name__ == "__main__":
    import sys
    
    in_file_name = sys.argv[1]
    plot_folder = sys.argv[2]
    
    main(in_file_name, plot_folder)

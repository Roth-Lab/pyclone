from __future__ import division

from clonal_estimation.independent_model import DataPoint, CollapsedLikelihood
from math import exp

import csv
import multiprocessing
import os

def main(in_file_name, plot_folder):
    reader = csv.DictReader(open(in_file_name), delimiter='\t')
    
    p = multiprocessing.Pool()
    
    p.map(run_gene, [x for x in reader])    
    
def run_gene(row):
    gene = row['gene']
    
    print gene
    
    a = int(row['ref_counts'])
    b = int(row['var_counts'])
    
    d = a + b
    
    cn = int(row['cn'])
    loh = int(row['loh'])
    
    mu, pi_r, pi_v = get_genotype_priors(cn, loh)
    
    data_point = DataPoint(a, d, pi_r, pi_v, mu)
    
    n = 1000
    likelihood = CollapsedLikelihood(data_point)
    
    x = []
    y = []
    
    norm_const = 0
    
    for i in range(1, n):
        x.append(i / n)
        y.append(exp(likelihood.compute_log_likelihood(x[-1])))
        
        norm_const += 1 / n * y[-1]
    
    y = [y_i / norm_const for y_i in y]
    
    out_file_name = os.path.join(plot_folder, '{0}.density.tsv'.format(gene))
    
    writer = csv.writer(open(out_file_name, 'w'), delimiter='\t')
    
    writer.writerow(x)
    writer.writerow(y) 

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

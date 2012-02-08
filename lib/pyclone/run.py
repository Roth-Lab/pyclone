'''
Created on 2012-02-08

@author: Andrew Roth
'''
from pyclone.model import DataPoint, BinomialLikelihood
from pyclone.samplers import DirichletProcessSampler

import cPickle
import csv

def run_dp_model(args):
    reader = csv.DictReader(open(args.in_file_name), delimiter='\t')

    genes = []
    likelihoods = []
    
    for row in reader:
        genes.append(row['gene'])
        
        a = int(row['a'])
        
        d = int(row['d'])
        
        mu_r = [float(x) for x in row['mu_r'].split(',')]
        delta_r = [float(x) for x in row['delta_r'].split(',')]
        
        mu_v = [float(x) for x in row['mu_v'].split(',')]
        delta_v = [float(x) for x in row['delta_v'].split(',')]
        
        data_point = DataPoint(a, d, mu_r, mu_v, delta_r, delta_v)

        likelihood = BinomialLikelihood(data_point)
        
        likelihoods.append(likelihood)
        
    db = {}
    
    db['input_file'] = args.in_file_name    
    db['genes'] = genes
    db['a'] = [x.a for x in likelihoods]
    db['d'] = [x.d for x in likelihoods]
    
    sampler = DirichletProcessSampler(likelihoods)
    
    sampler_results = sampler.sample(num_iters=args.num_iters, burnin=0, thin=1)
    
    db['results'] = sampler_results
    
    fh = open(args.out_file_name, 'wb')
    cPickle.dump(db, fh)
    fh.close()
'''
Created on 2012-02-08

@author: Andrew Roth
'''
from collections import OrderedDict

from pyclone.model import DataPoint, BinomialLikelihood
from pyclone.post_process import DpSamplerPostProcessor
from pyclone.results import AnalysisDB
from pyclone.samplers import DirichletProcessSampler

import csv
import os

def run_dp_model(args):
    '''
    Run a fresh instance of the DP model.
    '''
    data_set = load_data(args.in_file)
    
    analysis_db = AnalysisDB(args.out_prefix)
    
    analysis_db['input_file'] = open(args.in_file).readlines()
    
    analysis_db['genes'] = data_set.keys()
    
    analysis_db['data'] = data_set.values()
    
    sampler = DirichletProcessSampler(data_set.values(), burnin=args.burnin, thin=args.thin)
    
    sampler.sample(analysis_db, num_iters=args.num_iters)
    
    analysis_db['sampler'] = sampler
    
    analysis_db.close()

def resume_dp_model(args):
    '''
    Restart an existing analysis.
    '''    
    analysis_db = AnalysisDB(args.out_prefix)
    
    sampler = analysis_db['sampler']
    
    sampler.sample(analysis_db, num_iters=args.num_iters)
    
    analysis_db['sampler'] = sampler
    
    analysis_db.close()

def load_data(input_file_name):
    '''
    Load data from PyClone formatted input file.
    '''
    data = OrderedDict()
    
    reader = csv.DictReader(open(input_file_name), delimiter='\t')

    for row in reader:
        gene = row['gene']
        
        a = int(row['a'])
        
        d = int(row['d'])
        
        mu_r = [float(x) for x in row['mu_r'].split(',')]
        delta_r = [float(x) for x in row['delta_r'].split(',')]
        
        mu_v = [float(x) for x in row['mu_v'].split(',')]
        delta_v = [float(x) for x in row['delta_v'].split(',')]
        
        data[gene] = DataPoint(a, d, mu_r, mu_v, delta_r, delta_v)

    return data
        
def post_process_results(args):    
    analysis_db = AnalysisDB(args.results_prefix)
    
    post_processor = DpSamplerPostProcessor(analysis_db)
    
#    if args.raw_data:
#        raw_dir = os.path.join(args.out_dir, 'raw_data')
#        
#        safe_makedirs(raw_dir)
#        
#        write_raw_data(post_processor, raw_dir)
#    
#    if args.independent:
#        write_independent_posteriors(post_processor, args.out_dir)
    
    posteriors_dir = os.path.join(args.out_dir, 'posteriors')
    
    safe_makedirs(posteriors_dir)
    
    write_posteriors(post_processor, args.num_bins, posteriors_dir)
#
#def write_raw_data(post_processor, out_dir):
#    # Save genes
#    gene_file = os.path.join(out_dir, "genes.tsv")
#    writer = csv.writer(open(gene_file, 'w'), delimiter='\t')
#    writer.writerows(list_to_csv_rows(post_processor.genes))
#    
#    # Save alpha
#    alpha_file = os.path.join(out_dir, 'alpha.tsv')
#    writer = csv.writer(open(alpha_file, 'w'), delimiter='\t')
#    writer.writerows(list_to_csv_rows(post_processor.alpha))
#    
#    # Save num components
#    components_file = os.path.join(out_dir, 'components.tsv')
#    writer = csv.writer(open(components_file, 'w'), delimiter='\t')
#    writer.writerows(list_to_csv_rows(post_processor.num_components))
#    
#    # Save cellular frequencies.
#    cellular_freq_dir = os.path.join(out_dir, 'cellular_frequencies')
#    
#    if not os.path.exists(cellular_freq_dir):
#        os.makedirs(cellular_freq_dir)
#    
#    cellular_freqs = post_processor.cellular_frequencies
#    
#    for gene in post_processor.genes:
#        gene_file = os.path.join(cellular_freq_dir, "{0}.tsv".format(gene))
#        
#        fh = open(gene_file, 'w')
#        
#        writer = csv.writer(fh, delimiter='\t')
#        
#        writer.writerows(list_to_csv_rows(cellular_freqs[gene]))
#        
#        fh.close()
#    
#    # Save similarity matrix
#    sim_mat_file = os.path.join(out_dir, "similarity_matrix.tsv")
#    writer = csv.writer(open(sim_mat_file, 'w'), delimiter='\t')
#    writer.writerows(post_processor.similarity_matrix)
#
#def write_independent_posteriors(post_processor, out_dir):
#    pass
#
def write_posteriors(post_processor, num_bins, out_dir):
    # Save genes
    gene_file = os.path.join(out_dir, "genes.tsv")
    writer = csv.writer(open(gene_file, 'w'), delimiter='\t')
    writer.writerows(list_to_csv_rows(post_processor.genes))
    
    # Save alpha
    alpha_file = os.path.join(out_dir, 'alpha.tsv')
    writer = csv.writer(open(alpha_file, 'w'), delimiter='\t')
    writer.writerows(histogram_to_csv_rows(post_processor.get_alpha_posteriors()))
    
    # Save num components
    components_file = os.path.join(out_dir, 'components.tsv')
    writer = csv.writer(open(components_file, 'w'), delimiter='\t')
    writer.writerows(histogram_to_csv_rows(post_processor.get_num_component_posteriors()))
    
    # Save cellular frequencies.
    cellular_freq_dir = os.path.join(out_dir, 'cellular_frequencies')
    
    if not os.path.exists(cellular_freq_dir):
        os.makedirs(cellular_freq_dir)
    
    cellular_freqs = post_processor.get_cellular_frequency_posteriors(num_bins)
    
    for gene in post_processor.genes:
        gene_file = os.path.join(cellular_freq_dir, "{0}.tsv".format(gene))
        
        fh = open(gene_file, 'w')
        
        writer = csv.writer(fh, delimiter='\t')
        
        writer.writerows(histogram_to_csv_rows(cellular_freqs[gene]))
        
        fh.close()
    
    # Save similarity matrix
    sim_mat_file = os.path.join(out_dir, "similarity_matrix.tsv")
    writer = csv.writer(open(sim_mat_file, 'w'), delimiter='\t')
    writer.writerows(post_processor.get_similarity_posteriors())

def list_to_csv_rows(x):
    return [[x_i, ] for x_i in x]

def histogram_to_csv_rows(x):
    return zip(x[0], x[1])

def safe_makedirs(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)    
    

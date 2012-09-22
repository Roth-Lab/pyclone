'''
Created on 2012-02-08

@author: Andrew Roth
'''
from math import lgamma as log_gamma

from pyclone.sampler import DirichletProcessSampler, PyCloneData
from pyclone.trace import TraceDB, TracePostProcessor

import csv
import os

def run_dp_model(args):
    '''
    Run a fresh instance of the DP model.
    '''
    data, genes = load_pyclone_data(args.in_file)
    
    tace_db = TraceDB(args.trace_file, mode='w')
    
    tace_db['input_file'] = open(args.in_file).readlines()
    
    tace_db['genes'] = genes
    
    tace_db['data'] = data  
    
    sampler = DirichletProcessSampler(data, alpha=args.concentration)
    
    sampler.sample(tace_db, num_iters=args.num_iters)
    
    tace_db['sampler'] = sampler
    
    tace_db.close()

def resume_dp_model(args):
    '''
    Restart an existing analysis.
    '''    
    trace_db = TraceDB(args.trace_file, 'a')
    
    sampler = trace_db['sampler']
    
    sampler.sample(trace_db, num_iters=args.num_iters)
    
    trace_db['sampler'] = sampler
    
    trace_db.close()

def load_pyclone_data(file_name):
    '''
    Load data from PyClone formatted input file.
    '''
    data = []
    genes = []
    
    reader = csv.DictReader(open(file_name), delimiter='\t')

    for row in reader:
        genes.append(row['gene'])
        
        a = int(row['a'])
        
        d = int(row['d'])
        
        mu_r = [float(x) for x in row['mu_r'].split(',')]
        mu_v = [float(x) for x in row['mu_v'].split(',')]
        
        delta_r = [float(x) for x in row['delta_r'].split(',')]
        delta_v = [float(x) for x in row['delta_v'].split(',')]
        
        log_pi_r = get_log_mix_weights(delta_r)
        log_pi_v = get_log_mix_weights(delta_v)
        
        data.append(PyCloneData(a, d, tuple(mu_r), tuple(mu_v), tuple(log_pi_r), tuple(log_pi_v)))

    return data, genes

def get_log_mix_weights(delta):
    log_denominator = log_gamma(sum(delta) + 1)
    
    log_mix_weights = []
    
    for i, d_i in enumerate(delta):
        log_numerator = log_gamma(d_i + 1)
        
        for j, d_j in enumerate(delta):
            if i != j:
                log_numerator += log_gamma(d_j)
        
        log_mix_weights.append(log_numerator - log_denominator)
    
    return log_mix_weights
        
def post_process_trace(args):    
    safe_makedirs(args.out_dir)
    
    post_processor = TracePostProcessor(args.trace_file, args.burnin, args.thin)

    _write_trace(post_processor, args.out_dir)

def _write_trace(post_processor, out_dir):
    # Save alpha
    alpha_file = os.path.join(out_dir, 'alpha.tsv')
    writer = csv.writer(open(alpha_file, 'w'), delimiter='\t')
    writer.writerows(list_to_csv_rows(post_processor.alpha))
    
    # Save num components
    components_file = os.path.join(out_dir, 'num_components.tsv')
    writer = csv.writer(open(components_file, 'w'), delimiter='\t')
    writer.writerows(list_to_csv_rows(post_processor.num_components))
    
    # Save cellular frequencies.
    cellular_frequencies_file = os.path.join(out_dir, 'cellular_frequencies.tsv')    
    cellular_frequencies = post_processor.cellular_frequencies
    _write_trace_dict(cellular_frequencies, cellular_frequencies_file)
    
    # Save labels
    labels_file = os.path.join(out_dir, "labels.tsv")
    labels = post_processor.labels
    _write_trace_dict(labels, labels_file) 

def _write_trace_dict(trace_dict, file_name):
    fields = ['gene', 'trace']
    
    writer = csv.DictWriter(open(file_name, 'w'), fields, delimiter='\t')
    
    writer.writeheader()
    
    for gene in trace_dict:
        out_row = {
                   'gene' : gene,
                   'trace' : list_to_string(trace_dict[gene])
                   }
        
        writer.writerow(out_row)

def list_to_csv_rows(x):
    return [[x_i, ] for x_i in x]

def list_to_string(x):
    return ",".join([str(x_i) for x_i in x])

def histogram_to_csv_rows(x):
    return zip(x[0], x[1])

def safe_makedirs(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)    
    

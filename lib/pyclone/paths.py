'''
Created on Nov 28, 2015

@author: Andrew Roth
'''
import os
import yaml


def load_config(file_name):
    with open(file_name) as fh:
        config = yaml.load(fh)

    return config


def get_error_rates(config_file):
    values = {}

    config = load_config(config_file)

    for sample_id in get_sample_ids(config_file):
        values[sample_id] = config['samples'][sample_id]['error_rate']

    return values


def get_mutations_files(config_file):
    files = {}

    config = load_config(config_file)

    for sample_id in get_sample_ids(config_file):
        files[sample_id] = config['samples'][sample_id]['mutations_file']

    return files


def get_sample_ids(config_file):
    config = load_config(config_file)

    return config['samples'].keys()


def get_tumour_contents(config_file):
    values = {}

    config = load_config(config_file)

    for sample_id in get_sample_ids(config_file):
        values[sample_id] = config['samples'][sample_id]['tumour_content']['value']

    return values


def get_cellular_prevalence_trace_files(config_file):
    trace_files = {}

    trace_dir = get_trace_dir(config_file)

    for sample_id in get_sample_ids(config_file):
        trace_files[sample_id] = os.path.join(trace_dir, '{0}.cellular_prevalence.tsv.bz2'.format(sample_id))

    return trace_files


def get_concentration_trace_file(config_file):
    trace_dir = get_trace_dir(config_file)

    return os.path.join(trace_dir, 'alpha.tsv.bz2')


def get_labels_trace_file(config_file):
    trace_dir = get_trace_dir(config_file)

    return os.path.join(trace_dir, 'labels.tsv.bz2')


def get_precision_trace_file(config_file):
    trace_dir = get_trace_dir(config_file)

    return os.path.join(trace_dir, 'precision.tsv.bz2')


def get_trace_dir(config_file):
    config = load_config(config_file)

    return os.path.join(config['working_dir'], config['trace_dir'])

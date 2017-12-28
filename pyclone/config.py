'''
Created on 2013-02-12

@author: Andrew Roth
'''
from collections import namedtuple, OrderedDict
from math import log

import numpy as np
import yaml

try:
    from yaml import CLoader as Loader

except ImportError:
    from yaml import Loader

import pyclone.paths as paths

#=======================================================================================================================
# Load data for sampler
#=======================================================================================================================
PyCloneData = namedtuple(
    'PyCloneData',
    [
        'b',
        'd',
        'tumour_content',
        'cn_n',
        'cn_r',
        'cn_v',
        'mu_n',
        'mu_r',
        'mu_v',
        'log_pi'
    ]
)


class PyCloneConfig(object):
    def __init__(self, file_name):
        self._config = load_yaml_file(file_name)

        self._init_mutations()

        self._init_data()

        if len(self.mutations) == 0:
            raise Exception(' '.join((
                ('No mutations found in common across samples.'),
                ('This commonly occurs when the muation_id field does not match for mutations in the input files.'),
                ('Search the user group before posting a message or a bug.'),
            )))

    @property
    def base_measure_params(self):
        """ The parameters for the DP Beta base measure.
        """
        return self._config['base_measure_params']

    @property
    def beta_binomial_precision_prior(self):
        """ Precision parameter for the Beta-Binomial distribution.
        """
        return self._config['beta_binomial_precision'].get('prior', None)

    @property
    def beta_binomial_precision_proposal_precision(self):
        """ Precision of Beta proposal distribution for Metropolis-Hastings updates of Beta-Binomial precision.
        """
        return self._config['beta_binomial_precision'].get('proposal', {'precision': 0.1})['precision']

    @property
    def beta_binomial_precision_value(self):
        """ The (initial) value of the precision parameter for the Beta-Binomial distribution.
        """
        return self._config['beta_binomial_precision']['value']

    @property
    def concetration_prior(self):
        """ The prior for the DP concentration parameter.

        Return None if not prior specified.
        """
        return self._config['concentration'].get('prior', None)

    @property
    def concetration_value(self):
        """ The (initial) value of the DP concentration parameter.
        """
        return self._config['concentration']['value']

    @property
    def data(self):
        """ Data for all samples.
        """
        return self._data

    @property
    def density(self):
        return self._config['density']

    @property
    def error_rates(self):
        """ Sequencing error rates for each sample.
        """
        result = {}

        for sample_id in self.samples:
            result[sample_id] = self._config['samples'][sample_id]['error_rate']

        return result

    @property
    def init_method(self):
        """ Initialisation method for DP sampler.
        """
        return self._config.get('init_method', 'connected')

    @property
    def mutations(self):
        """ List of mutations in data.
        """
        return self._mutations

    @property
    def mutations_files(self):
        """ Dictionary of mutations files for samples.
        """
        result = {}

        for sample in self.samples:
            result[sample] = self._config['samples'][sample]['mutations_file']

        return result

    @property
    def samples(self):
        """ List of samples.
        """
        return list(self._config['samples'].keys())

    @property
    def tumour_content_values(self):
        """ Dictionary of tumour contents for each sample.
        """
        result = {}

        for sample_id in self.samples:
            result[sample_id] = self._config['samples'][sample_id]['tumour_content']['value']

        return result

    def _init_data(self):
        sample_data = OrderedDict()

        for sample_id in self.samples:
            sample_data[sample_id] = self._load_sample_data(sample_id)

        data = OrderedDict()

        for mutation_id in self.mutations:
            data[mutation_id] = OrderedDict()

            for sample_id in self.samples:
                data[mutation_id][sample_id] = sample_data[sample_id][mutation_id]

        self._data = data

    def _init_mutations(self):
        mutations = []

        for file_name in self.mutations_files.values():
            sample_config = load_yaml_file(file_name)

            sample_mutations = set()

            for mutation_config in sample_config['mutations']:
                sample_mutations.add(mutation_config['id'])

            mutations.append(sample_mutations)

        mutations = set.intersection(*mutations)

        self._mutations = sorted(mutations)

    def _load_sample_data(self, sample):
        """ Load data from PyClone YAML formatted input file.
        """
        data = OrderedDict()

        config = paths.load_config(self.mutations_files[sample])

        for mutation_dict in config['mutations']:
            mutation = load_mutation_from_dict(mutation_dict)

            data[mutation.id] = _get_pyclone_data(
                mutation, self.error_rates[sample], self.tumour_content_values[sample]
            )

        return data


def _get_pyclone_data(mutation, error_rate, tumour_content):
    a = mutation.ref_counts
    b = mutation.var_counts

    d = a + b

    cn_n = np.array([x.cn_n for x in mutation.states])
    cn_r = np.array([x.cn_r for x in mutation.states])
    cn_v = np.array([x.cn_v for x in mutation.states])

    mu_n = np.array([x.get_mu_n(error_rate) for x in mutation.states])
    mu_r = np.array([x.get_mu_r(error_rate) for x in mutation.states])
    mu_v = np.array([x.get_mu_v(error_rate) for x in mutation.states])

    prior_weights = tuple([x.prior_weight for x in mutation.states])

    log_pi = _get_log_pi(prior_weights)

    return PyCloneData(b, d, tumour_content, cn_n, cn_r, cn_v, mu_n, mu_r, mu_v, log_pi)


def _get_log_pi(weights):
    pi = [x / sum(weights) for x in weights]

    return np.array([log(x) for x in pi])

#=======================================================================================================================
# Parse mutation dict
#=======================================================================================================================


def get_mutation(mutation_id, ref_counts, var_counts, normal_cn, minor_cn, major_cn, prior):
    if minor_cn > minor_cn:
        raise Exception('{} is invalid. Major CN must be greater than minor CN.'.format(mutation_id))

    if major_cn == 0:
        raise Exception('{} is invalid. Major CN must be greater 0.'.format(mutation_id))

    states = _get_states(normal_cn, minor_cn, major_cn, prior)

    mutation = Mutation(mutation_id, ref_counts, var_counts)

    for (g_n, g_r, g_v) in states:
        state = State(g_n, g_r, g_v, 1)

        mutation.add_state(state)

    return mutation


def _get_states(normal_cn, minor_cn, major_cn, prior):
    if prior == 'major':
        states = _get_major_copy_states(normal_cn, minor_cn, major_cn)

    elif prior == 'parental':
        states = _get_parental_copy_number_states(normal_cn, minor_cn, major_cn)

    elif prior == 'total':
        states = _get_total_copy_number_states(normal_cn, minor_cn, major_cn)

    else:
        raise Exception('{0} is not a recognised method for setting variant priors.'.format(prior))

    return states


def _get_major_copy_states(normal_cn, minor_cn, major_cn):
    states = set()

    g_n = 'A' * normal_cn

    total_cn = minor_cn + major_cn

    for b in range(1, major_cn + 1):
        a = total_cn - b

        g_r = g_n

        g_v = 'A' * a + 'B' * b

        states.add((g_n, g_r, g_v))

    states.add((g_n, 'A' * total_cn, 'A' * (total_cn - 1) + 'B'))

    return sorted(states)


def _get_parental_copy_number_states(normal_cn, minor_cn, major_cn):
    states = set()

    g_n = 'A' * normal_cn

    total_cn = minor_cn + major_cn

    if normal_cn == total_cn:
        if minor_cn == 0:
            ref_genotypes = [g_n, g_n]

            var_genotypes = [
                'A' * minor_cn + 'B' * major_cn,
                'A' * (total_cn - 1) + 'B'
            ]
        else:
            ref_genotypes = [g_n, g_n, g_n]

            var_genotypes = [
                'A' * minor_cn + 'B' * major_cn,
                'A' * major_cn + 'B' * minor_cn,
                'A' * (total_cn - 1) + 'B'
            ]

    else:
        if minor_cn == 0:
            ref_genotypes = [g_n, g_n, 'A' * total_cn]

            var_genotypes = [
                'A' * minor_cn + 'B' * major_cn,
                'A' * (total_cn - 1) + 'B',
                'A' * (total_cn - 1) + 'B'
            ]
        else:
            ref_genotypes = [g_n, g_n, g_n, 'A' * total_cn]

            var_genotypes = [
                'A' * minor_cn + 'B' * major_cn,
                'A' * major_cn + 'B' * minor_cn,
                'A' * (total_cn - 1) + 'B',
                'A' * (total_cn - 1) + 'B'
            ]

    for g_r, g_v in zip(ref_genotypes, var_genotypes):
        state = (g_n, g_r, g_v)

        states.add(state)

    return sorted(states)


def _get_total_copy_number_states(normal_cn, minor_cn, major_cn):
    total_cn = minor_cn + major_cn

    states = set()

    g_n = 'A' * normal_cn

    var_genotypes = []

    for num_var_alleles in range(1, total_cn + 1):
        g_v = 'A' * (total_cn - num_var_alleles) + 'B' * num_var_alleles

        var_genotypes.append(g_v)

    ref_genotypes = set([g_n, 'A' * total_cn])

    for g_r in ref_genotypes:
        for g_v in var_genotypes:
            state = (g_n, g_r, g_v)

            states.add(state)

    return sorted(states)

#=======================================================================================================================
# Helper classes
#=======================================================================================================================


class Mutation(object):

    def __init__(self, mutation_id, ref_counts, var_counts):
        self.id = mutation_id

        self.ref_counts = ref_counts

        self.var_counts = var_counts

        self.states = []

    @property
    def cn_n(self):
        return [x.cn_n for x in self.states]

    @property
    def cn_r(self):
        return [x.cn_r for x in self.states]

    @property
    def cn_v(self):
        return [x.cn_v for x in self.states]

    @property
    def prior_weights(self):
        return [x.prior_weight for x in self.states]

    def add_state(self, state):
        self.states.append(state)

    def get_mu_n(self, error_rate):
        return [x.get_mu_n(error_rate) for x in self.states]

    def get_mu_r(self, error_rate):
        return [x.get_mu_r(error_rate) for x in self.states]

    def get_mu_v(self, error_rate):
        return [x.get_mu_v(error_rate) for x in self.states]

    def to_dict(self):
        return {
            'id': self.id,
            'ref_counts': self.ref_counts,
            'var_counts': self.var_counts,
            'states': [x.to_dict() for x in self.states]
        }


class State(object):

    def __init__(self, g_n, g_r, g_v, prior_weight):
        self.g_n = g_n

        self.g_r = g_r

        self.g_v = g_v

        self.prior_weight = prior_weight

    @property
    def cn_n(self):
        return len(self.g_n)

    @property
    def cn_r(self):
        return len(self.g_r)

    @property
    def cn_v(self):
        return len(self.g_v)

    def get_mu_n(self, error_rate):
        return self._get_variant_allele_probability(self.g_n, error_rate)

    def get_mu_r(self, error_rate):
        return self._get_variant_allele_probability(self.g_r, error_rate)

    def get_mu_v(self, error_rate):
        return self._get_variant_allele_probability(self.g_v, error_rate)

    def to_dict(self):
        return {'g_n': self.g_n, 'g_r': self.g_r, 'g_v': self.g_v, 'prior_weight': self.prior_weight}

    def _get_copy_number(self, genotype):
        if genotype is None:
            return 0
        else:
            return len(genotype)

    def _get_variant_allele_probability(self, genotype, error_rate):
        if genotype is None:
            return error_rate

        num_ref_alleles = genotype.count("A")
        num_var_alleles = genotype.count("B")

        cn = len(genotype)

        if cn != num_ref_alleles + num_var_alleles:
            raise Exception("{0} is not a valid genotype. Only A or B are allowed as alleles.")

        if num_ref_alleles == 0:
            return 1 - error_rate
        elif num_var_alleles == 0:
            return error_rate
        else:
            return num_var_alleles / cn

#=======================================================================================================================
# Factory functions
#=======================================================================================================================


def load_mutation_from_dict(d):
    mutation_id = d['id']

    ref_counts = int(d['ref_counts'])
    var_counts = int(d['var_counts'])

    mutation = Mutation(mutation_id, ref_counts, var_counts)

    for state_dict in d['states']:
        state = load_state_from_dict(state_dict)

        mutation.add_state(state)

    return mutation


def load_state_from_dict(d):
    g_n = d['g_n']
    g_r = d['g_r']
    g_v = d['g_v']

    prior_weight = float(d['prior_weight'])

    return State(g_n, g_r, g_v, prior_weight)


def load_yaml_file(file_name):
    with open(file_name) as fh:
        config = yaml.load(fh, Loader=Loader)

    return config

'''
Created on 2013-02-12

@author: Andrew Roth
'''
from collections import OrderedDict

import numpy as np
import yaml

try:
    from yaml import CLoader as Loader

except ImportError:
    from yaml import Loader

import pyclone.data
import pyclone.math_utils


class PyCloneConfig(object):
    def __init__(
            self,
            data_df,
            density='beta-binomial',
            grid_size=None,
            init_concentration=1.0,
            init_precision=400,
            over_ride_file=None,
            update_concentration=True,
            update_precision=True):

        self._config = {
            'base_measure_params': {'a': 1.0, 'b': 1.0},
            'beta_binomial_precision': {
                'prior': {'rate': 0.001, 'shape': 1.0},
                'proposal': {'precision': 0.01},
                'value': init_precision,
                'update': update_precision
            },
            'concentration': {
                'prior': {'rate': 0.1, 'shape': 0.1},
                'value': init_concentration,
                'update': update_concentration
            },
            'density': density,
            'grid_size': grid_size,
            'init_method': 'connected'
        }

        if over_ride_file is not None:
            with open(over_ride_file, 'r') as fh:
                self._config.update(yaml.load(fh))

        self._init_data(data_df)

    @staticmethod
    def from_trace(trace):
        new = PyCloneConfig.__new__(PyCloneConfig)

        new._config = trace.config

        new._init_data(trace['data'])

        return new

    @property
    def base_measure_params(self):
        """ The parameters for the DP Beta base measure.
        """
        return self._config['base_measure_params']

    @property
    def beta_binomial_precision_prior(self):
        """ Precision parameter for the Beta-Binomial distribution.
        """
        return self._config['beta_binomial_precision']['prior']

    @property
    def beta_binomial_precision_proposal_precision(self):
        """ Precision of Beta proposal distribution for Metropolis-Hastings updates of Beta-Binomial precision.
        """
        return self._config['beta_binomial_precision']['proposal']['precision']

    @property
    def beta_binomial_precision_value(self):
        """ The (initial) value of the precision parameter for the Beta-Binomial distribution.
        """
        return self._config['beta_binomial_precision']['value']

    @property
    def concentration_prior(self):
        """ The prior for the DP concentration parameter.

        Return None if not prior specified.
        """
        return self._config['concentration']['prior']

    @property
    def concentration_value(self):
        """ The (initial) value of the DP concentration parameter.
        """
        return self._config['concentration']['value']

    @property
    def discrete_approximation(self):
        return self._config['grid_size'] is not None

    @property
    def density(self):
        return self._config['density']

    @property
    def grid_size(self):
        return self._config['grid_size']

    @property
    def init_method(self):
        """ Initialisation method for DP sampler.
        """
        return self._config['init_method']

    @property
    def update_concentration(self):
        return self._concentration['concentration']['update']

    @property
    def update_precision(self):
        return self._config['beta_binomial_precision']['update']

    def to_dict(self):
        return self._config.copy()

    def update(self, config):
        self._config.update(config)

    def _init_data(self, df):
        self.samples = sorted(df['sample'].unique())

        # Filter for mutations present in all samples
        df = df.groupby(by='mutation').filter(lambda x: sorted(x['sample'].unique()) == self.samples)

        self.mutations = df['mutation'].unique()

        if len(self.mutations) == 0:
            raise Exception(' '.join((
                ('No mutations found in common across samples.'),
                ('This commonly occurs when the muation_id field does not match for mutations in the input files.'),
                ('Search the user group before posting a message or a bug.'),
            )))

        else:
            print('Samples: {}'.format(' '.join(self.samples)))

            print('Num mutations: {}'.format(len(self.mutations)))

        if 'error_rate' not in df.columns:
            df.loc[:, 'error_rate'] = 1e-3

        if 'tumour_content' not in df.columns:
            print('Tumour content column not found. Setting values to 1.0.')

            df.loc[:, 'tumour_content'] = 1.0

        self.data = OrderedDict()

        for name, mut_df in df.groupby('mutation'):
            sample_data_points = []

            mut_df = mut_df.set_index('sample')

            for sample in self.samples:
                row = mut_df.loc[sample]

                a = row['ref_counts']

                b = row['alt_counts']

                cn, mu, log_pi = get_major_cn_prior(
                    row['major_cn'],
                    row['minor_cn'],
                    row['normal_cn'],
                    error_rate=row['error_rate']
                )

                sample_data_points.append(pyclone.data.SampleDataPoint(a, b, cn, mu, log_pi, row['tumour_content']))

            self.data[name] = pyclone.data.DataPoint(self.samples, sample_data_points)


def get_major_cn_prior(major_cn, minor_cn, normal_cn, error_rate=1e-3):
    total_cn = major_cn + minor_cn

    cn = []

    mu = []

    log_pi = []

    # Consider all possible mutational genotypes consistent with mutation before CN change
    for x in range(1, major_cn + 1):
        cn.append((normal_cn, normal_cn, total_cn))

        mu.append((error_rate, error_rate, min(1 - error_rate, x / total_cn)))

        log_pi.append(0)

    # Consider mutational genotype of mutation before CN change if not already added
    mutation_after_cn = (normal_cn, total_cn, total_cn)

    if mutation_after_cn not in cn:
        cn.append(mutation_after_cn)

        mu.append((error_rate, error_rate, min(1 - error_rate, 1 / total_cn)))

        log_pi.append(0)

        assert len(set(cn)) == 2

    cn = np.array(cn, dtype=int)

    mu = np.array(mu, dtype=float)

    log_pi = pyclone.math_utils.log_normalize(np.array(log_pi, dtype=float))

    return cn, mu, log_pi

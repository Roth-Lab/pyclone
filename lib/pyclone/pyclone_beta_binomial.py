'''
Created on 2013-04-28

@author: Andrew Roth
'''
from __future__ import division

from collections import OrderedDict

from pydp.base_measures import BetaBaseMeasure, GammaBaseMeasure
from pydp.data import GammaData
from pydp.densities import Density
from pydp.proposal_functions import GammaProposal
from pydp.samplers.atom import BaseMeasureAtomSampler
from pydp.samplers.dp import DirichletProcessSampler
from pydp.samplers.global_params import MetropolisHastingsGlobalParameterSampler
from pydp.samplers.partition import AuxillaryParameterPartitionSampler

import numpy as np

from pyclone.math_utils import log_beta_binomial_likelihood, log_sum_exp, jit
from pyclone.multi_sample import MultiSampleBaseMeasure, MultiSampleDensity, MultiSampleAtomSampler
from pyclone.trace import DiskTrace

import pyclone.config as config


def run_pyclone_beta_binomial_analysis(config_file, num_iters, alpha, alpha_priors):
    data, sample_ids = config.load_data(config_file)

    sample_atom_samplers = OrderedDict()

    sample_base_measures = OrderedDict()

    sample_cluster_densities = OrderedDict()

    base_measure_params = config.load_base_measure_params(config_file)

    precision_params = config.load_precision_params(config_file)

    init_method = config.load_init_method(config_file)

    for sample_id in sample_ids:
        sample_base_measures[sample_id] = BetaBaseMeasure(base_measure_params['alpha'], base_measure_params['beta'])

        sample_cluster_densities[sample_id] = PyCloneBetaBinomialDensity(GammaData(precision_params['value']))

        sample_atom_samplers[sample_id] = BaseMeasureAtomSampler(
            sample_base_measures[sample_id],
            sample_cluster_densities[sample_id],
        )

    base_measure = MultiSampleBaseMeasure(sample_base_measures)

    cluster_density = MultiSampleDensity(sample_cluster_densities, shared_params=True)

    atom_sampler = MultiSampleAtomSampler(base_measure, cluster_density, sample_atom_samplers)

    partition_sampler = AuxillaryParameterPartitionSampler(base_measure, cluster_density)

    if 'prior' in precision_params:
        global_params_sampler = MetropolisHastingsGlobalParameterSampler(
            GammaBaseMeasure(
                precision_params['prior']['shape'],
                precision_params['prior']['rate']
            ),
            cluster_density,
            GammaProposal(precision_params['proposal']['precision'])
        )

    else:
        global_params_sampler = None

    sampler = DirichletProcessSampler(
        atom_sampler,
        partition_sampler,
        alpha,
        alpha_priors,
        global_params_sampler,
    )

    trace = DiskTrace(config_file, data.keys(), {'cellular_frequencies': 'x'}, precision=True)

    trace.open()

    sampler.sample(data.values(), trace, num_iters, init_method=init_method)

    trace.close()


class PyCloneBetaBinomialDensity(Density):

    def log_p(self, data, params):
        return self._log_p(data, params)

    def _log_p(self, data, params):
        return _log_p(
            data.b, data.d,
            data.cn_n, data.cn_r, data.cn_v,
            data.mu_n, data.mu_r, data.mu_v,
            data.log_pi,
            params.x, data.tumour_content, self.params.x
        )


@jit(cache=True, nopython=True)
def _log_p(b, d, cn_n, cn_r, cn_v, mu_n, mu_r, mu_v, log_pi, f, t, s):
    num_states = len(log_pi)

    ll = np.zeros(num_states)

    for i in range(num_states):
        ll[i] = log_pi[i]

        ll[i] += _log_beta_binomial_likelihood(
            b, d,
            cn_n[i], cn_r[i], cn_v[i],
            mu_n[i], mu_r[i], mu_v[i],
            f, t, s
        )

    return log_sum_exp(ll)


@jit(cache=True, nopython=True)
def _log_beta_binomial_likelihood(b, d, cn_n, cn_r, cn_v, mu_n, mu_r, mu_v, f, t, precision):
    p_n = (1 - t) * cn_n
    p_r = t * (1 - f) * cn_r
    p_v = t * f * cn_v

    norm_const = p_n + p_r + p_v

    p_n = p_n / norm_const
    p_r = p_r / norm_const
    p_v = p_v / norm_const

    mu = p_n * mu_n + p_r * mu_r + p_v * mu_v

    param_a = mu * precision

    param_b = (1 - mu) * precision

    return log_beta_binomial_likelihood(b, d, param_a, param_b)

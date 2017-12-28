'''
Created on 2013-04-28

@author: Andrew Roth
'''
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


def load_sampler(config):
    sample_atom_samplers = OrderedDict()

    sample_base_measures = OrderedDict()

    sample_cluster_densities = OrderedDict()

    for sample_id in config.samples:
        sample_base_measures[sample_id] = BetaBaseMeasure(**config.base_measure_params)

        sample_cluster_densities[sample_id] = PyCloneBetaBinomialDensity(
            GammaData(config.beta_binomial_precision_value)
        )

        sample_atom_samplers[sample_id] = BaseMeasureAtomSampler(
            sample_base_measures[sample_id],
            sample_cluster_densities[sample_id],
        )

    base_measure = MultiSampleBaseMeasure(sample_base_measures)

    cluster_density = MultiSampleDensity(sample_cluster_densities, shared_params=True)

    atom_sampler = MultiSampleAtomSampler(base_measure, cluster_density, sample_atom_samplers)

    partition_sampler = AuxillaryParameterPartitionSampler(base_measure, cluster_density)

    if config.beta_binomial_precision_prior is not None:
        global_params_sampler = MetropolisHastingsGlobalParameterSampler(
            GammaBaseMeasure(**config.beta_binomial_precision_prior),
            cluster_density,
            GammaProposal(config.beta_binomial_precision_proposal_precision)
        )

    else:
        global_params_sampler = None

    sampler = DirichletProcessSampler(
        atom_sampler,
        partition_sampler,
        config.concetration_value,
        config.concetration_prior,
        global_params_sampler,
    )

    return sampler


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

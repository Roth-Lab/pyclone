from collections import OrderedDict

from pydp.base_measures import BetaBaseMeasure, GammaBaseMeasure
from pydp.proposal_functions import GammaProposal
from pydp.samplers.atom import BaseMeasureAtomSampler
from pydp.samplers.dp import DirichletProcessSampler
from pydp.samplers.global_params import MetropolisHastingsGlobalParameterSampler
from pydp.samplers.partition import AuxillaryParameterPartitionSampler

import numpy as np

import pyclone.densities
import pyclone.multi_sample


def get_sampler(config):
    sample_atom_samplers = OrderedDict()

    sample_base_measures = OrderedDict()

    sample_cluster_densities = OrderedDict()

    for sample_id in config.samples:
        sample_base_measures[sample_id] = BetaBaseMeasure(**config.base_measure_params)

        if config.density == 'binomial':
            sample_cluster_densities[sample_id] = pyclone.densities.PyCloneBinomialDensity()

        else:
            sample_cluster_densities[sample_id] = pyclone.densities.PyCloneBetaBinomialDensity(
                config.beta_binomial_precision_value
            )

        sample_atom_samplers[sample_id] = BaseMeasureAtomSampler(
            sample_base_measures[sample_id],
            sample_cluster_densities[sample_id],
        )

    base_measure = pyclone.multi_sample.MultiSampleBaseMeasure(sample_base_measures)

    cluster_density = pyclone.multi_sample.MultiSampleDensity(sample_cluster_densities, shared_params=True)

    atom_sampler = pyclone.multi_sample.MultiSampleAtomSampler(base_measure, cluster_density, sample_atom_samplers)

    partition_sampler = AuxillaryParameterPartitionSampler(base_measure, cluster_density)

    if (config.density == 'binomial') or (config.beta_binomial_precision_prior is None):
        global_params_sampler = None

    else:
        global_params_sampler = MetropolisHastingsGlobalParameterSampler(
            GammaBaseMeasure(**config.beta_binomial_precision_prior),
            cluster_density,
            GammaProposal(config.beta_binomial_precision_proposal_precision)
        )

    sampler = DirichletProcessSampler(
        atom_sampler,
        partition_sampler,
        config.concetration_value,
        config.concetration_prior,
        global_params_sampler,
    )

    sampler.initialise_partition(config.init_method, len(config.data))

    return sampler


def run_mcmc(config, num_iters, sampler, trace):
    print('Beginning analysis using:')
    print('{} mutations'.format(len(config.mutations)))
    print('{} sample(s)'.format(len(config.samples)))
    print()

    for i in range(num_iters):
        sampler.interactive_sample(list(config.data.values()))

        state = sampler.state

        trace.update(state)

        if i % 100 == 0:
            print('Iteration: {}'.format(i))
            print('Number of clusters: {}'.format(len(np.unique(state['labels']))))
            print('DP concentration: {}'.format(state['alpha']))
            if state['global_params'] is not None:
                print('Beta-Binomial precision: {}'.format(state['global_params'][0]))
            print()

from collections import OrderedDict

from pgsm.mcmc.collapsed_gibbs import CollapsedGibbsSampler
from pgsm.mcmc.concentration import GammaPriorConcentrationSampler
from pgsm.mcmc.particle_gibbs_split_merge import ParticleGibbsSplitMergeSampler
from pgsm.partition_priors import DirichletProcessPartitionPrior
from pgsm.mcmc.split_merge_setup import UniformSplitMergeSetupKernel

from pydp.base_measures import BetaBaseMeasure, GammaBaseMeasure
from pydp.proposal_functions import GammaProposal
from pydp.samplers.atom import BaseMeasureAtomSampler
from pydp.samplers.dp import DirichletProcessSampler
from pydp.samplers.global_params import MetropolisHastingsGlobalParameterSampler
from pydp.samplers.partition import AuxillaryParameterPartitionSampler

import pgsm.distributions.pyclone
import numpy as np

import pyclone.densities
import pyclone.multi_sample


def get_sampler(config):
    if config.discrete_approximation:
        sampler = _load_marginal_sampler(config)

    else:
        sampler = _load_classic_sampler(config)

    return sampler


def _load_classic_sampler(config):
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
        config.concentration_value,
        config.concentration_prior,
        global_params_sampler,
    )

    sampler.initialise_partition(config.init_method, len(config.data))

    return sampler


def _load_marginal_sampler(config):
    return MarginalSampler(
        list(config.data.values()),
        concetration_prior=config.concentration_prior,
        init_concentration=config.concentration_value
    )


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
            if state.get('global_params', None) is not None:
                print('Beta-Binomial precision: {}'.format(state['global_params'][0]))
            print()


class MarginalSampler(object):
    def __init__(self, data, concetration_prior={'rate': 1.0, 'shape': 1.0}, init_concentration=1.0):
        self.dist = pgsm.distributions.pyclone.PyCloneDistribution(data[0].shape)

        self.partition_prior = DirichletProcessPartitionPrior(init_concentration)

        self.gibbs_sampler = CollapsedGibbsSampler(self.dist, self.partition_prior)

        self.setup_kernel = UniformSplitMergeSetupKernel(data, self.dist, self.partition_prior)

        self.pgsm_sampler = ParticleGibbsSplitMergeSampler.create_from_dist(
            self.dist, self.partition_prior, self.setup_kernel, num_anchors=2
        )

        self.conc_sampler = GammaPriorConcentrationSampler(concetration_prior['shape'], concetration_prior['rate'])

        self.pred_clustering = np.zeros(len(data))

    @property
    def state(self):
        return {
            'alpha': self.partition_prior.alpha,
            'labels': self.pred_clustering,
        }

    @state.setter
    def state(self, value):
        self.partition_prior.alpha = value['alpha']

        self.pred_clustering = value['labels']

    def interactive_sample(self, data):
        data = np.array(data)

        self.pred_clustering = self.pgsm_sampler.sample(self.pred_clustering, data)

        self.pred_clustering = self.gibbs_sampler.sample(self.pred_clustering, data)

        num_clusters = len(np.unique(self.pred_clustering))

        num_data_points = len(data)

        self.partition_prior.alpha = self.conc_sampler.sample(self.partition_prior.alpha, num_clusters, num_data_points)

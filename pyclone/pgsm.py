from collections import defaultdict, OrderedDict
from scipy.special import logsumexp as log_sum_exp

from pgsm.math_utils import discrete_rvs
from pgsm.mcmc.collapsed_gibbs import CollapsedGibbsSampler
from pgsm.mcmc.concentration import GammaPriorConcentrationSampler
from pgsm.mcmc.particle_gibbs_split_merge import ParticleGibbsSplitMergeSampler
from pgsm.mcmc.split_merge_setup import UniformSplitMergeSetupKernel
from pgsm.partition_priors import DirichletProcessPartitionPrior

import numpy as np
import pandas as pd

import pyclone.math_utils


class MarginalSampler(object):
    def __init__(self, config):
        self.config = config

        self._init_data(config)

        self._sampler = self._init_sampler(config)

    @property
    def ccfs(self):
        pass

    @property
    def state(self):
        return {
            'alpha': self.partition_prior.alpha,
            'labels': self.pred_clustering,
            'params': self._sample_data_params()
        }

    @state.setter
    def state(self, value):
        self.partition_prior.alpha = value['alpha']

        self.pred_clustering = value['labels']

    def interactive_sample(self):
        self.pred_clustering = self.pgsm_sampler.sample(self.pred_clustering, self.data)

        self.pred_clustering = self.gibbs_sampler.sample(self.pred_clustering, self.data)

        num_clusters = len(np.unique(self.pred_clustering))

        num_data_points = len(self.data)

        if self.conc_sampler is not None:
            self.partition_prior.alpha = self.conc_sampler.sample(
                self.partition_prior.alpha, num_clusters, num_data_points
            )

    def _init_data(self, config):
        data = []

        for data_point in config.data.values():
            data.append(
                data_point.to_likelihood_grid(
                    config.density,
                    config.grid_size,
                    precision=config.beta_binomial_precision_value
                )
            )

        self.data = np.array(data)

    def _init_sampler(self, config):
        self.dist = PyCloneDistribution(self.data[0].shape)

        self.partition_prior = DirichletProcessPartitionPrior(config.concentration_value)

        self.gibbs_sampler = CollapsedGibbsSampler(self.dist, self.partition_prior)

        self.setup_kernel = UniformSplitMergeSetupKernel(self.data, self.dist, self.partition_prior)

        self.pgsm_sampler = ParticleGibbsSplitMergeSampler.create_from_dist(
            self.dist, self.partition_prior, self.setup_kernel, num_anchors=2
        )

        if config.concentration_prior is None:
            self.conc_sampler = None

        else:
            self.conc_sampler = GammaPriorConcentrationSampler(
                config.concentration_prior['shape'], config.concentration_prior['rate']
            )

        if config.init_method == 'connected':
            self.pred_clustering = np.zeros(len(self.data))

        elif config.init_method == 'disconnected':
            self.pred_clustering = np.arange(len(self.data))

        else:
            raise Exception('Unknown initialisation method {}.'.format(config.init_method))

    def _sample_cluster_params(self):
        result = defaultdict(OrderedDict)

        ccfs = np.linspace(0, 1, self.data[0].shape[1])

        dist = PyCloneDistribution(self.data[0].shape)

        for cluster_idx in np.unique(self.pred_clustering):
            cluster_data = self.data[self.pred_clustering == cluster_idx]

            cluster_params = dist.create_params_from_data(cluster_data)

            log_posterior = cluster_params.normalized_log_pdf_grid

            for s_idx, sample in enumerate(self.config.samples):
                ccf_idx = discrete_rvs(log_posterior[s_idx, :])

                result[cluster_idx][sample] = ccfs[ccf_idx]

        return result

    def _sample_data_params(self):
        cluster_params = self._sample_cluster_params()

        params = OrderedDict()

        for sample in self.config.samples:
            params[sample] = [cluster_params[cluster_idx][sample] for cluster_idx in self.pred_clustering]

        params = pd.DataFrame(params, index=self.config.mutations)

        return params


class PycloneParameters(object):

    def __init__(self, log_pdf_grid, N):
        self.log_pdf_grid = log_pdf_grid

        self.N = N

    @property
    def normalized_log_pdf_grid(self):
        result = []

        for i in range(self.log_pdf_grid.shape[0]):
            result.append(pyclone.math_utils.log_normalize(self.log_pdf_grid[i, :]))

        return np.array(result)

    def copy(self):
        return PycloneParameters(self.log_pdf_grid.copy(), self.N)

    def decrement(self, x):
        self.log_pdf_grid -= x

        self.N -= 1

    def increment(self, x):
        self.log_pdf_grid += x

        self.N += 1


class PyCloneDistribution(object):

    def __init__(self, grid_size):
        self.grid_size = grid_size

    def create_params(self):
        return PycloneParameters(np.zeros(self.grid_size), 0)

    def create_params_from_data(self, X):
        X = np.atleast_3d(X)

        return PycloneParameters(np.sum(X, axis=0), X.shape[0])

    def log_marginal_likelihood(self, params):
        log_p = np.sum(
            log_sum_exp(params.log_pdf_grid - np.log(self.grid_size[1]), axis=1)
        )

        return log_p

    def log_predictive_likelihood(self, data_point, params):
        params.increment(data_point)

        ll = self.log_marginal_likelihood(params)

        params.decrement(data_point)

        ll -= self.log_marginal_likelihood(params)

        return ll

    def log_predictive_likelihood_bulk(self, data, params):
        log_p = np.zeros(len(data))

        for i, data_point in enumerate(data):
            log_p[i] = self.log_predictive_likelihood(data_point, params)

        return log_p

    def log_pairwise_marginals(self, data, params):
        num_data_points = len(data)

        log_p = np.zeros((num_data_points, num_data_points))

        for i in range(num_data_points):
            for j in range(num_data_points):
                if i == j:
                    continue

                params = self.create_params_from_data(data[[i, j]])

                log_p[i, j] = self.log_marginal_likelihood(params)

        return log_p

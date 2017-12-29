from scipy.special import logsumexp as log_sum_exp
from pgsm.mcmc.collapsed_gibbs import CollapsedGibbsSampler
from pgsm.mcmc.concentration import GammaPriorConcentrationSampler
from pgsm.mcmc.particle_gibbs_split_merge import ParticleGibbsSplitMergeSampler
from pgsm.partition_priors import DirichletProcessPartitionPrior
from pgsm.mcmc.split_merge_setup import UniformSplitMergeSetupKernel

import numpy as np

import pyclone.math_utils


def load_marginal_sampler(config):
    return MarginalSampler(
        list(config.data.values()),
        concetration_prior=config.concentration_prior,
        init_concentration=config.concentration_value
    )


class MarginalSampler(object):
    def __init__(self, data, concetration_prior={'rate': 1.0, 'shape': 1.0}, init_concentration=1.0):
        self.dist = PyCloneDistribution(data[0].shape)

        self.partition_prior = DirichletProcessPartitionPrior(init_concentration)

        self.gibbs_sampler = CollapsedGibbsSampler(self.dist, self.partition_prior)

        self.setup_kernel = UniformSplitMergeSetupKernel(data, self.dist, self.partition_prior)

        self.pgsm_sampler = ParticleGibbsSplitMergeSampler.create_from_dist(
            self.dist, self.partition_prior, self.setup_kernel, num_anchors=2
        )

        if concetration_prior is None:
            self.conc_sampler = None

        else:
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

        if self.conc_sampler is not None:
            self.partition_prior.alpha = self.conc_sampler.sample(
                self.partition_prior.alpha, num_clusters, num_data_points)


class PycloneParameters(object):

    def __init__(self, log_pdf_grid, N):
        self.log_pdf_grid = log_pdf_grid

        self.N = N

    @property
    def normalized_log_pdf_grid(self):
        return pyclone.math_utils.log_normalize(self.log_pdf_grid)

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

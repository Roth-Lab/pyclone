from scipy.special import logsumexp as log_sum_exp

import numpy as np

import pyclone.math_utils


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
        log_p = np.sum(log_sum_exp(params, axis=1))
        
        print(log_sum_exp(params, axis=1))

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

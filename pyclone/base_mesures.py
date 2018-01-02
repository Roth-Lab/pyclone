from pydp.base_measures import BaseMeasure, BetaBaseMeasure
from pydp.data import BetaData
from scipy.misc import logsumexp as log_sum_exp

import numpy as np
import scipy.stats
from pgsm.math_utils import discrete_rvs


class PyCloneBaseMeasure(object):
    def as_pgsm(self, grid_size):
        raise NotImplementedError

    def as_pydp(self):
        raise NotImplementedError


class BetaPyCloneBaseMeasure(PyCloneBaseMeasure):
    def __init__(self, a, b):
        self.a = a

        self.b = b

    def as_pgsm(self, grid_size):
        return scipy.stats.beta.logpdf(np.linspace(0, 1, grid_size))

    def as_pydp(self):
        return BetaBaseMeasure(self.a, self.b)


class PointMassPyCloneBaseMeasure(PyCloneBaseMeasure):
    def __init__(self, a, b, pi):
        self.a = a

        self.b = b

        self.log_pi = np.log(np.array(pi) / np.sum(pi))

    def as_pgsm(self, grid_size):
        p_0 = np.ones(grid_size) * float('-inf')

        p_0[0] = self.log_pi[0]

        p_1 = np.ones(grid_size) * float('-inf')

        p_1[-1] = self.log_pi[1]

        p_sc = self.log_pi[2] + scipy.stats.beta.logpdf(np.linspace(0, 1, grid_size), self.a, self.b)

        return log_sum_exp([p_0, p_1, p_sc], axis=0)

    def as_pydp(self):
        return PointMassBaseMeasure(self.a, self.b, self.log_pi)


class PointMassBaseMeasure(BaseMeasure):

    def __init__(self, a, b, log_pi):
        self.a = a

        self.b = b

        self.log_pi = log_pi

        pi = np.exp(log_pi)

        self.pi = pi / np.sum(pi)

    def log_p(self, data):
        log_p = self.log_pi[2] + scipy.stats.beta.logpdf(data.x, self.a, self.b)

        if data.x == 0:
            log_p += self.log_pi[0]

        elif data.x == 1:
            log_p += self.log_pi[1]

        return log_p

    def random(self):
        idx = discrete_rvs(self.pi)

        if idx == 0:
            x = 0.0

        elif idx == 1:
            x = 1.0

        else:
            x = scipy.stats.beta.rvs(self.a, self.b)

        return BetaData(x)

from collections import OrderedDict

import numba
import numpy as np

import pyclone.math_utils


class DataPoint(object):
    def __init__(self, samples, sample_data_points):
        self.samples = samples

        self.sample_data_points = sample_data_points

    def get_ccf_grid(self, grid_size):
        return np.linspace(0, 1, grid_size)

    def to_dict(self):
        return OrderedDict(zip(self.samples, self.sample_data_points))

    def to_likelihood_grid(self, density, grid_size, precision=None):
        if (density == 'beta-binomial') and (precision is None):
            raise Exception('Precision must be set when using Beta-Binomial.')

        shape = (len(self.samples), grid_size)

        log_ll = np.zeros(shape)

        for s_idx, data_point in enumerate(self.sample_data_points):
            for i, ccf in enumerate(self.get_ccf_grid(grid_size)):
                if density == 'beta-binomial':
                    log_ll[s_idx, i] = pyclone.math_utils.log_pyclone_beta_binomial_pdf(data_point, ccf, precision)

                elif density == 'binomial':
                    log_ll[s_idx, i] = pyclone.math_utils.log_pyclone_binomial_pdf(data_point, ccf)

        return log_ll


@numba.jitclass([
    ('a', numba.int64),
    ('b', numba.int64),
    ('cn', numba.int64[:, :]),
    ('mu', numba.float64[:, :]),
    ('log_pi', numba.float64[:]),
    ('t', numba.float64)
])
class SampleDataPoint(object):
    def __init__(self, a, b, cn, mu, log_pi, t):
        self.a = a
        self.b = b
        self.cn = cn
        self.mu = mu
        self.log_pi = log_pi
        self.t = t

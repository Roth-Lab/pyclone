from pydp.data import GammaData
from pydp.densities import Density

import pyclone.math_utils


class PyCloneBetaBinomialDensity(Density):
    def __init__(self, params):
        self.params = GammaData(params)

    def log_p(self, data, params):
        return self._log_p(data, params)

    def _log_p(self, data, params):
        return pyclone.math_utils.log_pyclone_beta_binomial_pdf(data, params.x, self.params.x)


class PyCloneBinomialDensity(Density):
    def log_p(self, data, params):
        return self._log_p(data, params)

    def _log_p(self, data, params):
        return pyclone.math_utils.log_pyclone_binomial_pdf(data, params.x)

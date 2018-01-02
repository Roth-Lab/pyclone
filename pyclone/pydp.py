from collections import OrderedDict, namedtuple

from pydp.base_measures import BaseMeasure, BetaBaseMeasure, GammaBaseMeasure
from pydp.data import BetaData, GammaData
from pydp.densities import Density
from pydp.partition import PartitionCell
from pydp.proposal_functions import GammaProposal
from pydp.proposal_functions import ProposalFunction
from pydp.samplers.atom import AtomSampler, BaseMeasureAtomSampler
from pydp.samplers.dp import DirichletProcessSampler
from pydp.samplers.global_params import MetropolisHastingsGlobalParameterSampler
from pydp.samplers.partition import AuxillaryParameterPartitionSampler

import pandas as pd

import pyclone.math_utils


class InstantiatedSampler(object):
    def __init__(self, config):
        self.config = config

        self._init_data(config)

        self._sampler = self._init_sampler(config)

    @property
    def state(self):
        state = self._sampler.state

        params = OrderedDict()

        for sample in self.config.samples:
            params[sample] = [data_point_param[sample].x for data_point_param in state['params']]

        state['params'] = pd.DataFrame(params, index=self.config.mutations)

        return state

    @state.setter
    def state(self, value):
        params = []

        for _, row in value['params'].iterrows():
            row = row.apply(lambda x: BetaData(x))

            params.append(row.to_dict())

        value['params'] = params

        self._sampler.state = value

    def interactive_sample(self):
        self._sampler.interactive_sample(self.data)

    def _init_data(self, config):
        self.data = []

        for data_point in config.data.values():
            self.data.append(data_point.to_dict())

    def _init_sampler(self, config):
        base_measure = pyclone.pydp.MultiSampleBaseMeasure.from_samples(
            BetaBaseMeasure, config.base_measure_params, config.samples
        )

        if config.density == 'binomial':
            density_cls = pyclone.pydp.PyCloneBinomialDensity

            density_params = {}

        elif config.density == 'beta-binomial':
            density_cls = pyclone.densities.PyCloneBetaBinomialDensity

            density_params = {'params': config.beta_binomial_precision_value}

        cluster_density = pyclone.pydp.MultiSampleDensity.from_samples(
            density_cls, density_params, config.samples
        )

        atom_sampler = pyclone.pydp.MultiSampleAtomSampler.from_samples(
            BaseMeasureAtomSampler, base_measure, cluster_density, config.samples
        )

        partition_sampler = AuxillaryParameterPartitionSampler(base_measure, cluster_density)

        if (config.density == 'binomial') or (config.update_precision is None):
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


class MultiSampleAtomSampler(AtomSampler):

    def __init__(self, base_measure, cluster_density, atom_samplers):
        AtomSampler.__init__(self, base_measure, cluster_density)

        self.atom_samplers = atom_samplers

    @staticmethod
    def from_samples(atom_sampler_cls, base_measure, cluster_density, samples):
        atom_samplers = OrderedDict()

        for sample in samples:
            atom_samplers[sample] = atom_sampler_cls(
                base_measure.base_measures[sample], cluster_density.cluster_densities[sample]
            )

        return MultiSampleAtomSampler(base_measure, cluster_density, atom_samplers)

    def sample_atom(self, data, cell):
        new_atom = OrderedDict()

        for sample_id in self.atom_samplers:
            sample_data = [x[sample_id] for x in data]

            sample_cell = PartitionCell(cell.value[sample_id])

            sample_cell._items = cell._items

            new_atom[sample_id] = self.atom_samplers[sample_id].sample_atom(sample_data, sample_cell)

        return new_atom


class MultiSampleBaseMeasure(BaseMeasure):

    def __init__(self, base_measures):
        '''
        Args:
            base_measures: (dict) Mapping of sample IDs to base measures.
        '''
        self.base_measures = base_measures

    @staticmethod
    def from_samples(base_measure_cls, base_measure_params, samples):
        base_measures = OrderedDict()

        for sample in samples:
            base_measures[sample] = base_measure_cls(**base_measure_params)

        return MultiSampleBaseMeasure(base_measures)

    def log_p(self, data):
        log_p = 0

        for sample_id in self.base_measures:
            log_p += self.base_measures[sample_id].log_p(data[sample_id])

        return log_p

    def random(self):
        random_sample = OrderedDict()

        for sample_id in self.base_measures:
            random_sample[sample_id] = self.base_measures[sample_id].random()

        return random_sample


class MultiSampleDensity(Density):
    '''
    Wraps a collection of univariate densities.
    '''

    def __init__(self, cluster_densities, shared_params=False):
        '''
        Args:
            cluster_densities: (dict) A collection of Density objects for each sample.
        '''
        self.cluster_densities = cluster_densities

        self.shared_params = shared_params

    @staticmethod
    def from_samples(density_cls, density_params, samples):
        densities = OrderedDict()

        for sample in samples:
            densities[sample] = density_cls(**density_params)

        return MultiSampleDensity(densities, shared_params=True)

    @property
    def params(self):
        if self.shared_params:
            for cluster_id in self.cluster_densities:
                return self.cluster_densities[cluster_id].params

        else:
            params = OrderedDict()

            for cluster_id in self.cluster_densities:
                params[cluster_id] = self.cluster_densities[cluster_id].params

            return params

    @params.setter
    def params(self, value):
        if self.shared_params:
            for cluster_id in self.cluster_densities:
                self.cluster_densities[cluster_id].params = value

        elif isinstance(value, namedtuple):
            for cluster_id in self.cluster_densities:
                self.cluster_densities[cluster_id].params = value[cluster_id]

        else:
            raise Exception('Cannot set object type {0} as a density parameter'.format(type(value)))

    def log_p(self, data, params):
        log_p = 0

        for sample_id in self.cluster_densities:
            density = self.cluster_densities[sample_id]

            log_p += density.log_p(data[sample_id], params[sample_id])

        return log_p


class MultiSampleProposalFunction(ProposalFunction):

    def __init__(self, proposal_funcs):
        '''
        Args:
            proposal_funcs: (dict) A collection of ProposalFunction for each sample.
        '''
        self.proposal_funcs = proposal_funcs

    def log_p(self, data, params):
        log_p = 0

        for sample_id in self.proposal_funcs:
            log_p += self.proposal_funcs[sample_id].log_p(data[sample_id], params[sample_id])

        return log_p

    def random(self, params):
        random_sample = OrderedDict()

        for sample_id in self.proposal_funcs:
            random_sample[sample_id] = self.proposal_funcs[sample_id].random(params[sample_id])

        return random_sample


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

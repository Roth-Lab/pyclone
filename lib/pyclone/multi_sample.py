'''
Functions to allow multiple samples to be analysed in PyDP framework.

Created on 2013-04-28

@author: Andrew Roth
'''
from collections import OrderedDict, namedtuple
from pydp.base_measures import BaseMeasure
from pydp.densities import Density
from pydp.proposal_functions import ProposalFunction
from pydp.samplers.atom import AtomSampler
from pydp.partition import PartitionCell


class MultiSampleAtomSampler(AtomSampler):

    def __init__(self, base_measure, cluster_density, atom_samplers):
        AtomSampler.__init__(self, base_measure, cluster_density)

        self.atom_samplers = atom_samplers

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

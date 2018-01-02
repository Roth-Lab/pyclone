'''
Created on 2013-04-23

@author: Andrew Roth
'''
from pydp.data import GammaData

import os
import pandas as pd
import tables
import warnings

warnings.filterwarnings('ignore', category=pd.io.pytables.PerformanceWarning)
warnings.filterwarnings('ignore', category=pd.core.common.SettingWithCopyWarning)
warnings.filterwarnings('ignore', category=tables.exceptions.NaturalNameWarning)


class DiskTrace(object):

    def __init__(self, file_name, mode='r'):
        if mode == 'w':
            if os.path.exists(file_name):
                raise Exception('Trace file {} exists. Remove the file or specify another location.'.format(file_name))

        self._store = pd.HDFStore(file_name, complevel=9, mode=mode)

        if mode == 'w':
            self._idx = 0

        else:
            self._idx = len(self['alpha'])

    def __getitem__(self, key):
        return self._store[key]

    def __setitem__(self, key, value):
        self._store[key] = value

    @property
    def cancer_cell_fractions(self):
        result = {}

        for sample in self.samples:
            key = 'params/{}'.format(sample)

            result[sample] = self[key].pivot(values='value', columns='mutation', index='idx')

        return result

    @property
    def config(self):
        return self['config'].iloc[0]

    @config.setter
    def config(self, value):
        self['config'] = pd.Series([value])

    @property
    def labels(self):
        return self['labels'].pivot(values='value', columns='mutation', index='idx')

    @property
    def mutations(self):
        return self['mutations']

    @mutations.setter
    def mutations(self, value):
        self['mutations'] = pd.Series(value)

    @property
    def samples(self):
        return self['samples']

    @samples.setter
    def samples(self, value):
        self['samples'] = pd.Series(value)

    @property
    def state(self):
        state = {
            'alpha': self['alpha'].iloc[-1],
            'labels': self['/state/labels'],
            'params': self['/state/params']
        }

        if '/beta_binomial_precision' in self._store.keys():
            state['global_params'] = GammaData(self['beta_binomial_precision'].iloc[-1])

        else:
            state['global_params'] = None

        return state

    def close(self):
        self._store.close()

    def update(self, state):
        self._store.append('alpha', pd.Series({self._idx: state['alpha']}))

        if state.get('global_params', None) is not None:
            precision = float(state['global_params'].x)

            self._store.append('beta_binomial_precision', pd.Series({self._idx: precision}))

        labels = pd.DataFrame({'mutation': self.mutations, 'value': state['labels']})

        labels['idx'] = self._idx

        self._store.append('labels', labels)

        self['/state/labels'] = pd.Series(state['labels'], index=self.mutations)

        for sample in self.samples:
            sample_params = state['params'][sample].reset_index()

            sample_params.columns = 'mutation', 'value'

            sample_params['value'] = sample_params['value'].astype(float)

            sample_params['idx'] = self._idx

            self._store.append('params/{}'.format(sample), sample_params)

        self['/state/params'] = state['params']

        self._idx += 1

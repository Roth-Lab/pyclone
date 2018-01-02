'''
Created on 2013-04-23

@author: Andrew Roth
'''
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
            key = 'ccfs/{}'.format(sample)

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
            'ccfs': self['/state/ccfs'],
            'labels': self['/state/labels'].values
        }

        if '/beta_binomial_precision' in self._store.keys():
            state['beta_binomial_precision'] = self['beta_binomial_precision'].iloc[-1]

        else:
            state['beta_binomial_precision'] = None

        return state

    def close(self):
        self._store.close()

    def update(self, state):
        self._store.append('alpha', pd.Series({self._idx: state['alpha']}))

        if state.get('beta_binomial_precision', None) is not None:
            self._store.append('beta_binomial_precision', pd.Series({self._idx: state['beta_binomial_precision']}))

        labels = pd.DataFrame({'mutation': self.mutations, 'value': state['labels']})

        labels['idx'] = self._idx

        self._store.append('labels', labels)

        self['/state/labels'] = pd.Series(state['labels'], index=self.mutations)

        for sample in self.samples:
            sample_ccfs = state['ccfs'][sample].reset_index()

            sample_ccfs.columns = 'mutation', 'value'

            sample_ccfs['value'] = sample_ccfs['value'].astype(float)

            sample_ccfs['idx'] = self._idx

            self._store.append('ccfs/{}'.format(sample), sample_ccfs)

        self['/state/ccfs'] = state['ccfs']

        self._idx += 1

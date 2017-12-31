import numpy as np

import pyclone.pydp
import pyclone.pgsm


def get_sampler(config):
    if config.discrete_approximation:
        sampler = pyclone.pgsm.MarginalSampler(config)

    else:
        sampler = pyclone.pydp.InstantiatedSampler(config)

    return sampler


def run_mcmc(config, num_iters, sampler, trace):
    print('Beginning analysis using:')
    print('{} mutations'.format(len(config.mutations)))
    print('{} sample(s)'.format(len(config.samples)))
    print()

    for i in range(1, num_iters + 1):
        sampler.interactive_sample()

        state = sampler.state

        trace.update(state)

        if i % 10 == 0:
            print('Iteration: {}'.format(i))
            print('Number of clusters: {}'.format(len(np.unique(state['labels']))))
            print('DP concentration: {}'.format(state['alpha']))
            if state.get('global_params', None) is not None:
                print('Beta-Binomial precision: {}'.format(state['global_params'][0]))
            print()

'''
Created on Nov 29, 2015

@author: Andrew Roth
'''
from pydp.cluster import cluster_with_mpear
from pydp.data import BetaData, GammaData
from pydp.utils import log_space_normalise

import numpy as np
import pandas as pd

from pyclone.pyclone_beta_binomial import PyCloneBetaBinomialDensity
from pyclone.pyclone_binomial import PyCloneBinomialDensity


def cluster_pyclone_trace(trace, burnin=0, thin=1, max_clusters=None):
    X = trace.labels.values[burnin::thin]

    labels = cluster_with_mpear(X, max_clusters=max_clusters)

    labels = pd.Series(labels, index=trace.mutations)

    labels = labels.reset_index()

    labels.columns = 'mutation_id', 'cluster_id'

    return labels


def load_summary_table(config, trace, burnin=0, grid_size=101, max_clusters=None, min_size=0, thin=1):
    df = load_table(
        config, trace, burnin=burnin, grid_size=grid_size, max_clusters=max_clusters, min_size=min_size, thin=thin
    )

    df = df.set_index(['sample_id', 'cluster_id', 'size'])

    x = df.columns.astype(float).values[np.newaxis, :]

    p = np.exp(df.values)

    m_1 = np.sum(x * p, axis=1)

    m_2 = np.sum(np.power(x, 2) * p, axis=1)

    var = m_2 - np.power(m_1, 2)

    std = np.sqrt(var)

    mean_df = pd.DataFrame(m_1, index=df.index)

    std_df = pd.DataFrame(std, index=df.index)

    out_df = pd.concat([mean_df, std_df], axis=1)

    out_df.columns = 'mean', 'std'

    out_df = out_df.reset_index()

    return out_df


def load_table(config, trace, burnin=0, grid_size=101, min_size=0, max_clusters=None, thin=1):
    if config.density == 'pyclone_beta_binomial':
        precision = trace['beta_binomial_precision']

        precision = precision.iloc[burnin::thin].mean()

        density = PyCloneBetaBinomialDensity(GammaData(precision))

    elif config.density == 'pyclone_binomial':
        density = PyCloneBinomialDensity()

    else:
        raise Exception('Only pyclone_binomial and pyclone_beta_binomial density are supported.')

    data = config.data

    labels = cluster_pyclone_trace(trace, burnin=burnin, max_clusters=max_clusters, thin=thin)

    labels = labels.set_index('mutation_id')

    posteriors = []

    for cluster_id, cluster_df in labels.groupby('cluster_id'):
        mutation_ids = list(cluster_df.index)

        cluster_data = [data[x] for x in mutation_ids]

        for sample_id in config.samples:
            cluster_sample_data = [x[sample_id] for x in cluster_data]

            cluster_sample_posterior = _compute_posterior(cluster_sample_data, density, grid_size)

            cluster_sample_posterior['sample_id'] = sample_id

            cluster_sample_posterior['cluster_id'] = cluster_id

            cluster_sample_posterior['size'] = len(mutation_ids)

            posteriors.append(cluster_sample_posterior)

    df = pd.DataFrame(posteriors)

    df = df.set_index(['sample_id', 'cluster_id', 'size'])

    df = df.reset_index()

    df = df[df['size'] >= min_size]

    return df


def _compute_posterior(data, density, grid_size):
    posterior = {}

    for cellular_prevalence in np.linspace(0, 1, grid_size):
        posterior[cellular_prevalence] = 0

        for data_point in data:
            posterior[cellular_prevalence] += density.log_p(data_point, BetaData(cellular_prevalence))

    posterior = dict(list(zip(list(posterior.keys()), log_space_normalise(list(posterior.values())))))

    return posterior

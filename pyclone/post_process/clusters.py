'''
Created on Nov 29, 2015

@author: Andrew Roth
'''
from pydp.cluster import cluster_with_mpear

import numpy as np
import pandas as pd

import pyclone.math_utils


def cluster_pyclone_trace(trace, burnin=0, thin=1, max_clusters=None):
    X = trace.labels.values[burnin::thin]

    labels = cluster_with_mpear(X, max_clusters=max_clusters)

    labels = pd.Series(labels, index=trace.mutations)

    labels = labels.reset_index()

    labels.columns = 'mutation', 'cluster'

    return labels


def load_summary_table(config, trace, burnin=0, grid_size=101, max_clusters=None, min_size=0, thin=1):
    df = load_table(
        config, trace, burnin=burnin, grid_size=grid_size, max_clusters=max_clusters, min_size=min_size, thin=thin
    )

    df = df.set_index(['sample', 'cluster', 'size'])

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


def load_table(config, trace, burnin=0, grid_size=101, max_clusters=None, min_size=0, thin=1):
    if config.density == 'beta-binomial':
        if config.update_precision:
            precision = trace['beta_binomial_precision']

            precision = precision.iloc[burnin::thin].mean()

        else:
            precision = config.beta_binomial_precision_value

    elif config.density == 'binomial':
        precision = None

    else:
        raise Exception('Only pyclone_binomial and pyclone_beta_binomial density are supported.')

    labels = cluster_pyclone_trace(trace, burnin=burnin, max_clusters=max_clusters, thin=thin)

    labels = labels.set_index('mutation')

    ccfs = list(config.data.values())[0].get_ccf_grid(grid_size)

    posteriors = []

    for cluster_id, cluster_df in labels.groupby('cluster'):
        mutation_ids = list(cluster_df.index)

        cluster_data = [
            config.data[x].to_likelihood_grid(config.density, grid_size, precision=precision) for x in mutation_ids
        ]

        for s_idx, sample in enumerate(config.samples):
            cluster_sample_data = np.array([x[s_idx] for x in cluster_data])

            cluster_sample_posterior = pyclone.math_utils.log_normalize(cluster_sample_data.sum(axis=0))

            cluster_sample_posterior = dict(
                zip(ccfs, cluster_sample_posterior)
            )

            cluster_sample_posterior['sample'] = sample

            cluster_sample_posterior['cluster'] = cluster_id

            cluster_sample_posterior['size'] = len(mutation_ids)

            posteriors.append(cluster_sample_posterior)

    df = pd.DataFrame(posteriors)

    df = df.set_index(['sample', 'cluster', 'size'])

    df = df.reset_index()

    df = df[df['size'] >= min_size]

    return df

'''
Created on Nov 29, 2015

@author: Andrew Roth
'''
from pydp.cluster import cluster_with_mpear
from pydp.data import BetaData, GammaData
from pydp.utils import log_space_normalise

import numpy as np
import pandas as pd

from pyclone.config import load_data
from pyclone.pyclone_beta_binomial import PyCloneBetaBinomialDensity
from pyclone.pyclone_binomial import PyCloneBinomialDensity

import pyclone.paths as paths
import pyclone.trace as trace


def cluster_pyclone_trace(config_file, burnin, thin):
    labels_trace = trace.load_cluster_labels_trace(
        paths.get_labels_trace_file(config_file),
        burnin,
        thin
    )

    X = labels_trace.values

    labels = cluster_with_mpear(X)

    labels = pd.Series(labels, index=labels_trace.columns)

    labels = labels.reset_index()

    labels.columns = 'mutation_id', 'cluster_id'

    return labels


def load_summary_table(config_file, burnin=0, mesh_size=101, min_size=0, thin=1):

    df = load_table(
        config_file,
        burnin=burnin,
        mesh_size=mesh_size,
        min_size=min_size,
        thin=thin,
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


def load_table(config_file, burnin=0, min_size=0, mesh_size=101, thin=1):
    config = paths.load_config(config_file)

    if config['density'] == 'pyclone_beta_binomial':
        precision_file = paths.get_precision_trace_file(config_file)

        precision = pd.read_csv(precision_file, header=None, compression='bz2', sep='\t', squeeze=True)

        precision = precision.iloc[burnin::thin].mean()

        density = PyCloneBetaBinomialDensity(GammaData(precision))

    elif config['density'] == 'pyclone_binomial':
        density = PyCloneBinomialDensity()

    else:
        raise Exception('Only pyclone_binomial and pyclone_beta_binomial density are supported.')

    data, sample_ids = load_data(config_file)

    labels = cluster_pyclone_trace(config_file, burnin, thin)

    labels = labels.set_index('mutation_id')

    posteriors = []

    for cluster_id, cluster_df in labels.groupby('cluster_id'):
        mutation_ids = list(cluster_df.index)

        cluster_data = [data[x] for x in mutation_ids]

        for sample_id in sample_ids:
            cluster_sample_data = [x[sample_id] for x in cluster_data]

            cluster_sample_posterior = _compute_posterior(cluster_sample_data, density, mesh_size)

            cluster_sample_posterior['sample_id'] = sample_id

            cluster_sample_posterior['cluster_id'] = cluster_id

            cluster_sample_posterior['size'] = len(mutation_ids)

            posteriors.append(cluster_sample_posterior)

    df = pd.DataFrame(posteriors)

    df = df.set_index(['sample_id', 'cluster_id', 'size'])

    df = df.reset_index()

    df = df[df['size'] >= min_size]

    return df


def _compute_posterior(data, density, mesh_size):
    posterior = {}

    for cellular_prevalence in np.linspace(0, 1, mesh_size):
        posterior[cellular_prevalence] = 0

        for data_point in data:
            posterior[cellular_prevalence] += density.log_p(data_point, BetaData(cellular_prevalence))

    posterior = dict(zip(posterior.keys(), log_space_normalise(posterior.values())))

    return posterior

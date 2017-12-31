'''
Created on Nov 27, 2015

@author: Andrew Roth
'''
import pandas as pd

import pyclone.post_process.clusters


def load_table(config, trace, burnin, thin, max_clusters=None, min_cluster_size=0, old_style=False):
    data = _load_variant_allele_frequencies(trace)

    labels = pyclone.post_process.clusters.cluster_pyclone_trace(
        trace, burnin=burnin, max_clusters=max_clusters, thin=thin
    )

    labels = labels.set_index('mutation')['cluster']

    cluster_sizes = labels.value_counts()

    used_clusters = cluster_sizes[cluster_sizes >= min_cluster_size].index

    labels = labels[labels.isin(used_clusters)]

    labels = labels.reset_index()

    data = pd.merge(data, labels, on='mutation', how='inner')

    cols = [
        'mutation', 'sample', 'cluster',
        'ref_counts', 'alt_counts',
        'major_cn', 'minor_cn', 'normal_cn',
        'ccf', 'ccf_std', 'vaf'
    ]

    ccf = _load_cancer_cell_fractions(trace, burnin, thin)

    data = pd.merge(data, ccf, how='inner', on=['mutation', 'sample'])

    data = data[cols]

    data = data.sort_values(by=['cluster', 'mutation', 'sample'])

    if old_style:
        data = _reformat_multi_sample_table(data)

    return data


def _reformat_multi_sample_table(df):
    df = df.rename(
        columns={'cluster': 'cluster_id', 'mutation': 'mutation_id', 'sample': 'sample_id'}
    )

    mean_df = df[['mutation_id', 'sample_id', 'ccf']]

    mean_df = mean_df.pivot(index='mutation_id', columns='sample_id', values='ccf')

    std_df = df[['mutation_id', 'sample_id', 'ccf_std']]

    std_df = std_df.pivot(index='mutation_id', columns='sample_id', values='ccf_std')

    std_df = std_df.rename(columns=lambda x: '{0}_std'.format(x))

    cluster_df = df[['mutation_id', 'cluster_id']]

    cluster_df = cluster_df.groupby('mutation_id')['cluster_id'].apply(lambda x: x.iloc[0])

    return pd.concat([mean_df, std_df, cluster_df], axis=1).reset_index()


def _load_variant_allele_frequencies(trace):
    data = trace['data']

    data['vaf'] = data['alt_counts'] / (data['ref_counts'] + data['alt_counts'])

    return data[['mutation', 'sample', 'ref_counts', 'alt_counts', 'major_cn', 'minor_cn', 'normal_cn', 'vaf']]


def _load_cancer_cell_fractions(trace, burnin, thin):
    data = []

    for sample_id, sample_data in trace.cancer_cell_fractions.items():
        sample_data = sample_data.iloc[burnin::thin]

        sample_data = _load_sample_cancer_cell_fractions(sample_data)

        sample_data['sample'] = sample_id

        data.append(sample_data)

    data = pd.concat(data, axis=0)

    # Filter for mutations in all samples
    data = data[data['mutation'].isin(trace.mutations)]

    return data


def _load_sample_cancer_cell_fractions(data):
    data = pd.concat([data.mean(), data.std()], axis=1)

    data.columns = 'ccf', 'ccf_std'

    data.index.name = 'mutation'

    data = data.reset_index()

    return data

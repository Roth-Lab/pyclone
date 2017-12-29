'''
Created on Nov 27, 2015

@author: Andrew Roth
'''
import pandas as pd

import pyclone.post_process.clusters


def load_table(config, trace, burnin, thin, max_clusters=None, min_cluster_size=0, old_style=False):
    ccf = _load_cancer_cell_fractions(trace, burnin, thin)

    vaf = _load_variant_allele_frequencies(config)

    data = pd.merge(vaf, ccf, how='inner', on=['mutation_id', 'sample_id'])

    labels = pyclone.post_process.clusters.cluster_pyclone_trace(
        trace, burnin=burnin, max_clusters=max_clusters, thin=thin
    )

    labels = labels.set_index('mutation_id')['cluster_id']

    cluster_sizes = labels.value_counts()

    used_clusters = cluster_sizes[cluster_sizes >= min_cluster_size].index

    labels = labels[labels.isin(used_clusters)]

    labels = labels.reset_index()

    data = pd.merge(data, labels, on='mutation_id', how='inner')

    cols = ['mutation_id', 'sample_id', 'cluster_id', 'ccf', 'ccf_std', 'vaf']

    data = data[cols]

    data = data.sort_values(by=['cluster_id', 'mutation_id', 'sample_id'])

    if old_style:
        data = _reformat_multi_sample_table(data)

    return data


def _reformat_multi_sample_table(df):
    mean_df = df[['mutation_id', 'sample_id', 'ccf']]

    mean_df = mean_df.pivot(index='mutation_id', columns='sample_id', values='ccf')

    std_df = df[['mutation_id', 'sample_id', 'ccf_std']]

    std_df = std_df.pivot(index='mutation_id', columns='sample_id', values='ccf_std')

    std_df = std_df.rename(columns=lambda x: '{0}_std'.format(x))

    cluster_df = df[['mutation_id', 'cluster_id']]

    cluster_df = cluster_df.groupby('mutation_id')['cluster_id'].apply(lambda x: x.iloc[0])

    return pd.concat([mean_df, std_df, cluster_df], axis=1).reset_index()


def _load_variant_allele_frequencies(config):
    data = []

    for sample in config.samples:
        sample_data = _load_sample_variant_allele_frequencies(config, sample)

        sample_data['sample_id'] = sample

        data.append(sample_data)

    data = pd.concat(data, axis=0)

    # Filter for mutations in all samples
    data = data[data['mutation_id'].isin(config.mutations)]

    return data


def _load_sample_variant_allele_frequencies(config, sample):
    '''
    Load data from PyClone formatted input file.
    '''
    data = {}

    for mutation_id in config.data:
        mutation = config.data[mutation_id][sample]

        try:
            data[mutation_id] = mutation.b / (mutation.a + mutation.b)

        except ZeroDivisionError:
            data[mutation_id] = pd.np.nan

    data = pd.Series(data, name='vaf')

    data.index.name = 'mutation_id'

    data = data.reset_index()

    return data


def _load_cancer_cell_fractions(trace, burnin, thin):
    data = []

    for sample_id, sample_data in trace.cancer_cell_fractions.items():
        sample_data = sample_data.iloc[burnin::thin]

        sample_data = _load_sample_cancer_cell_fractions(sample_data)

        sample_data['sample_id'] = sample_id

        data.append(sample_data)

    data = pd.concat(data, axis=0)

    # Filter for mutations in all samples
    data = data[data['mutation_id'].isin(trace.mutations)]

    return data


def _load_sample_cancer_cell_fractions(data):
    data = pd.concat([data.mean(), data.std()], axis=1)

    data.columns = 'ccf', 'ccf_std'

    data.index.name = 'mutation_id'

    data = data.reset_index()

    return data

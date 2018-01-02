'''
Created on 2012-02-08

@author: Andrew Roth
'''
import numpy as np
import pandas as pd
import random

import pyclone.config
import pyclone.mcmc
import pyclone.post_process.clusters
import pyclone.post_process.plot
import pyclone.trace


def build_table(trace_file, out_file, table_format, burnin=0, grid_size=101, max_clusters=100, thin=1):
    trace = pyclone.trace.DiskTrace(trace_file)

    config = _load_config_from_trace(trace)

    if table_format == 'cluster':
        df = pyclone.post_process.clusters.load_summary_table(
            config,
            trace,
            burnin=burnin,
            grid_size=grid_size,
            max_clusters=max_clusters,
            thin=thin,
        )

    elif table_format == 'loci':
        df = pyclone.post_process.loci.load_table(
            config,
            trace,
            burnin=burnin,
            max_clusters=max_clusters,
            old_style=False,
            thin=thin
        )

    elif table_format == 'old':
        df = pyclone.post_process.loci.load_table(
            config,
            trace,
            burnin=burnin,
            max_clusters=max_clusters,
            old_style=True,
            thin=thin
        )

    df.to_csv(out_file, index=False, sep='\t')

    trace.close()


def plot_clusters(
        trace_file,
        out_file,
        burnin=0,
        grid_size=101,
        plot_format='density',
        max_clusters=100,
        min_cluster_size=0,
        samples=[],
        thin=1):

    trace = pyclone.trace.DiskTrace(trace_file)

    config = _load_config_from_trace(trace)

    kwargs = {
        'burnin': burnin,
        'grid_size': grid_size,
        'max_clusters': max_clusters,
        'min_cluster_size': min_cluster_size,
        'samples': samples,
        'thin': thin
    }

    if plot_format == 'density':
        pyclone.post_process.plot.clusters.density_plot(config, trace, out_file, **kwargs)

    elif plot_format == 'line':
        pyclone.post_process.plot.clusters.parallel_coordinates_plot(
            config, trace, out_file, **kwargs
        )

    elif plot_format == 'scatter':
        pyclone.post_process.plot.clusters.scatter_plot(
            config, trace, out_file, **kwargs
        )

    trace.close()


def plot_loci(
        trace_file,
        out_file,
        plot_format,
        burnin=0,
        max_clusters=None,
        min_cluster_size=0,
        samples=None,
        thin=1):

    trace = pyclone.trace.DiskTrace(trace_file)

    config = _load_config_from_trace(trace)

    kwargs = {
        'burnin': burnin,
        'max_clusters': max_clusters,
        'min_cluster_size': min_cluster_size,
        'samples': samples,
        'thin': thin
    }

    if plot_format.startswith('ccf'):
        kwargs['value'] = 'ccf'

    elif plot_format.startswith('vaf'):
        kwargs['value'] = 'vaf'

    if plot_format == 'ccf-density':
        [kwargs.pop(x) for x in list(kwargs.keys()) if 'cluster' in x]

        kwargs.pop('value')

        pyclone.post_process.plot.loci.density_plot(
            trace,
            out_file,
            **kwargs
        )

    elif plot_format in ['ccf-line', 'vaf-line']:
        pyclone.post_process.plot.loci.parallel_coordinates_plot(
            config,
            trace,
            out_file,
            **kwargs
        )

    elif plot_format in ['ccf-scatter', 'vaf-scatter']:
        pyclone.post_process.plot.loci.scatter_plot(
            config,
            trace,
            out_file,
            **kwargs
        )

    elif plot_format == 'similarity-matrix':
        kwargs.pop('samples')

        pyclone.post_process.plot.loci.similarity_matrix_plot(
            trace,
            out_file,
            **kwargs
        )

    trace.close()


def resume_analysis(num_iters, trace_file):
    trace = pyclone.trace.DiskTrace(trace_file, mode='a')

    config = _load_config_from_trace(trace)

    sampler = pyclone.mcmc.get_sampler(config)

    sampler.state = trace.state

    pyclone.mcmc.run_mcmc(config, num_iters, sampler, trace)

    trace.close()


def run_analysis(
        in_file,
        trace_file,
        concentration_value=1.0,
        config_file=None,
        density='beta-binomial',
        grid_size=None,
        no_concentration_update=False,
        no_precision_update=False,
        num_iters=int(1e4),
        precision_value=1.0,
        seed=None):

    if seed is not None:
        np.random.seed(seed)

        random.seed(seed)

    trace = pyclone.trace.DiskTrace(trace_file, mode='w')

    data_df = pd.read_csv(in_file, sep='\t')

    if grid_size is not None:
        no_precision_update = True

    config = pyclone.config.PyCloneConfig(
        data_df,
        density=density,
        grid_size=grid_size,
        init_concentration=concentration_value,
        init_precision=precision_value,
        over_ride_file=config_file,
        update_concentration=(not no_concentration_update),
        update_precision=(not no_precision_update)
    )

    trace.config = config.to_dict()

    trace['data'] = data_df

    trace.mutations = config.mutations

    trace.samples = config.samples

    sampler = pyclone.mcmc.get_sampler(config)

    pyclone.mcmc.run_mcmc(config, num_iters, sampler, trace)

    trace.close()


def _load_config_from_trace(trace):
    config = pyclone.config.PyCloneConfig.from_trace(trace)

    return config

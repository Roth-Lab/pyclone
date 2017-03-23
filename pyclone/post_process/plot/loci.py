'''
Created on Nov 30, 2015

@author: Andrew Roth
'''
from scipy.cluster.hierarchy import average
from scipy.spatial.distance import pdist, squareform

import matplotlib.gridspec as gs
import matplotlib.pyplot as pp
import pandas as pd
import seaborn as sb

import pyclone.paths as paths
import pyclone.post_process as post_process
import pyclone.trace as trace

import defaults
import _scatter
import utils

#=======================================================================================================================
# Density
#=======================================================================================================================


def density_plot(
        config_file,
        plot_file,
        burnin=0,
        samples=None,
        thin=1):

    utils.setup_plot()

    df = _load_density_df(config_file, burnin, thin)

    if samples is None:
        samples = sorted(df['sample_id'].unique())

    else:
        df = df[df['sample_id'].isin(samples)]

    loci = df['mutation_id'].unique()

    num_loci = len(loci)

    width = 8

    height = 2 * num_loci + 2

    fig = pp.figure(figsize=(width, height))

    grid = gs.GridSpec(nrows=num_loci, ncols=1)

    for ax_index, locus in enumerate(loci):
        ax = fig.add_subplot(grid[ax_index])

        utils.setup_axes(ax)

        plot_df = df[df['mutation_id'] == locus]

        sb.violinplot(
            ax=ax,
            data=plot_df,
            x='sample_id',
            y='cellular_prevalence',

            inner=None,
            order=samples,
            scale='width'
        )

        ax.set_ylabel('')

        if ax_index != (num_loci - 1):
            ax.set_xticklabels([])

            ax.set_xlabel('')

        else:
            ax.set_xlabel(defaults.sample_label)

        ax.set_ylim(*defaults.cellular_prevalence_limits)

        ax.annotate(
            locus,
            xy=(1.01, 0.5),
            xycoords='axes fraction',
            fontsize=defaults.axis_label_font_size
        )

        utils.set_tick_label_font_sizes(ax, defaults.tick_label_font_size)

    fig.text(
        -0.01,
        0.5,
        defaults.cellular_prevalence_label,
        fontsize=defaults.axis_label_font_size,
        ha='center',
        rotation=90,
        va='center'
    )

    grid.tight_layout(fig, h_pad=3)

    utils.save_figure(fig, plot_file)


def _load_density_df(config_file, burnin, thin):
    trace_files = paths.get_cellular_prevalence_trace_files(config_file)

    df = []

    for sample_id, file_name in trace_files.items():

        sample_df = trace.load_cellular_frequencies_trace(
            file_name,
            burnin,
            thin
        )

        sample_df = sample_df.unstack().reset_index()

        sample_df.columns = 'mutation_id', 'iter', 'cellular_prevalence'

        sample_df['sample_id'] = sample_id

        sample_df = sample_df[['mutation_id', 'sample_id', 'cellular_prevalence']]

        df.append(sample_df)

    df = pd.concat(df)

    return df

#=======================================================================================================================
# Parallel coordinates
#=======================================================================================================================


def parallel_coordinates_plot(
        config_file,
        plot_file,
        burnin=0,
        max_clusters=None,
        min_cluster_size=0,
        samples=None,
        thin=1,
        value='cellular_prevalence'):

    utils.setup_plot()

    df = post_process.loci.load_table(
        config_file,
        burnin,
        thin,
        max_clusters=max_clusters,
        min_cluster_size=min_cluster_size
    )

    color_map = utils.get_clusters_color_map(df['cluster_id'])

    if samples is None:
        samples = sorted(df['sample_id'].unique())

    else:
        df = df[df['sample_id'].isin(samples)]

    df['sample_index'] = df['sample_id'].apply(lambda x: samples.index(x))

    df = df.sort_values(by='sample_index')

    fig = pp.figure()

    ax = fig.add_subplot(1, 1, 1)

    utils.setup_axes(ax)

    for cluster_id, cluster_df in df.groupby('cluster_id'):
        for _, locus_df in cluster_df.groupby('mutation_id'):
            x = locus_df['sample_index']

            y = locus_df[value]

            ax.plot(
                x,
                y,
                alpha=0.75,
                c=color_map[cluster_id],
                marker=defaults.line_plot_marker,
                markersize=defaults.line_plot_marker_size
            )

    ax.set_xlabel(defaults.sample_label, fontsize=defaults.axis_label_font_size)

    if value == 'cellular_prevalence':
        ax.set_ylabel(defaults.cellular_prevalence_label, fontsize=defaults.axis_label_font_size)

    elif value == 'variant_allele_frequency':
        ax.set_ylabel(defaults.variant_allele_frequency_label)

    ax.set_xticks(sorted(df['sample_index'].unique()))

    ax.set_xticklabels(samples)

    utils.set_tick_label_font_sizes(ax, defaults.tick_label_font_size)

    ax.set_ylim(*defaults.cellular_prevalence_limits)

    box = ax.get_position()

    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    legend_handles = utils.get_legend_handles(color_map)

    legend = ax.legend(
        legend_handles.values(),
        legend_handles.keys(),
        bbox_to_anchor=(1, 0.5),
        fontsize=defaults.legend_font_size,
        loc='center left',
        title=defaults.cluster_label
    )

    legend.get_title().set_fontsize(defaults.legend_title_font_size)

    utils.save_figure(fig, plot_file)

#=======================================================================================================================
# Scatter
#=======================================================================================================================


def scatter_plot(
        config_file,
        plot_file,
        burnin=0,
        max_clusters=None,
        min_cluster_size=0,
        samples=None,
        thin=1,
        value='cellular_prevalence'):

    utils.setup_plot()

    df = post_process.loci.load_table(
        config_file,
        burnin,
        thin,
        max_clusters=max_clusters,
        min_cluster_size=min_cluster_size
    )

    if samples is None:
        samples = sorted(df['sample_id'].unique())

    color_map = utils.get_clusters_color_map(df['cluster_id'])

    cluster_df = df[['mutation_id', 'cluster_id']].drop_duplicates().set_index('mutation_id')

    loci_color_map = cluster_df['cluster_id'].map(color_map).to_dict()

    mean_df = df.pivot(index='mutation_id', columns='sample_id', values=value)

    _scatter.plot_all_pairs(
        loci_color_map,
        mean_df,
        plot_file,
        samples,
        legend_color_map=color_map
    )

#=======================================================================================================================
# Similarity matrix
#=======================================================================================================================


def similarity_matrix_plot(
        config_file,
        plot_file,
        burnin=0,
        max_clusters=None,
        min_cluster_size=0,
        samples=None,
        thin=1):

    sb.set_style('whitegrid')

    labels = post_process.cluster_pyclone_trace(config_file, burnin, thin, max_clusters=max_clusters)

    labels = labels.set_index('mutation_id')

    labels = labels['cluster_id']

    color_map = utils.get_clusters_color_map(labels)

    cluster_sizes = labels.value_counts()

    used_clusters = cluster_sizes[cluster_sizes >= min_cluster_size].index

    labels = labels[labels.isin(used_clusters)]

    used_loci = labels.index

    trace_file = paths.get_labels_trace_file(config_file)

    labels_trace = trace.load_cluster_labels_trace(trace_file, burnin, thin)

    labels_trace = labels_trace[used_loci]

    dist_mat = pdist(labels_trace.values.T, 'hamming')

    Z = average(dist_mat)

    dist_mat = pd.DataFrame(squareform(dist_mat), index=labels_trace.columns, columns=labels_trace.columns)

    sim_mat = 1 - dist_mat

    N = sim_mat.shape[0]

    cluster_colors = labels.map(color_map)

    size = 0.12 * N

    g = sb.clustermap(
        sim_mat,
        cmap='Blues',
        col_colors=cluster_colors,
        row_colors=cluster_colors,
        col_linkage=Z,
        row_linkage=Z,
        figsize=(size, size)
    )

    ax = g.ax_heatmap

    utils.set_tick_label_font_sizes(ax, defaults.small_tick_label_font_size)

    utils.set_tick_label_rotations(ax)

    ax.set_xlabel('Loci', fontsize=defaults.axis_label_font_size)

    ax.set_ylabel('Loci', fontsize=defaults.axis_label_font_size)

    g.fig.savefig(plot_file, bbox_inches='tight')

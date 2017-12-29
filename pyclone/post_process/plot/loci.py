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

import pyclone.post_process as post_process

from pyclone.post_process.plot import defaults
from pyclone.post_process.plot import scatter
from pyclone.post_process.plot import utils

#=======================================================================================================================
# Density
#=======================================================================================================================


def density_plot(
        trace,
        out_file,
        burnin=0,
        samples=None,
        thin=1):

    utils.setup_plot()

    df = _load_density_df(trace, burnin, thin)

    if len(samples) == 0:
        samples = sorted(df['sample'].unique())

    else:
        df = df[df['sample'].isin(samples)]

    loci = df['mutation'].unique()

    num_loci = len(loci)

    width = 8

    height = 2 * num_loci + 2

    fig = pp.figure(figsize=(width, height))

    grid = gs.GridSpec(nrows=num_loci, ncols=1)

    for ax_index, locus in enumerate(loci):
        ax = fig.add_subplot(grid[ax_index])

        utils.setup_axes(ax)

        plot_df = df[df['mutation'] == locus]

        sb.violinplot(
            ax=ax,
            data=plot_df,
            x='sample',
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

    utils.save_figure(fig, out_file)


def _load_density_df(trace, burnin=0, thin=1):
    df = []

    for sample, sample_df in trace.cancer_cell_fractions.items():
        sample_df = sample_df[burnin::thin]

        sample_df = sample_df.unstack().reset_index()

        sample_df.columns = 'mutation', 'iter', 'cellular_prevalence'

        sample_df['sample'] = sample

        sample_df = sample_df[['mutation', 'sample', 'cellular_prevalence']]

        df.append(sample_df)

    df = pd.concat(df)

    return df

#=======================================================================================================================
# Parallel coordinates
#=======================================================================================================================


def parallel_coordinates_plot(
        config,
        trace,
        out_file,
        burnin=0,
        max_clusters=None,
        min_cluster_size=0,
        samples=None,
        thin=1,
        value='ccf'):

    utils.setup_plot()

    df = post_process.loci.load_table(
        config,
        trace,
        burnin=burnin,
        thin=thin,
        max_clusters=max_clusters,
        min_cluster_size=min_cluster_size
    )

    color_map = utils.get_clusters_color_map(df['cluster'])

    if len(samples) == 0:
        samples = sorted(df['sample'].unique())

    else:
        df = df[df['sample'].isin(samples)]

    df['sample_index'] = df['sample'].apply(lambda x: samples.index(x))

    df = df.sort_values(by='sample_index')

    fig = pp.figure()

    ax = fig.add_subplot(1, 1, 1)

    utils.setup_axes(ax)

    for cluster, cluster_df in df.groupby('cluster'):
        for _, locus_df in cluster_df.groupby('mutation'):
            x = locus_df['sample_index']

            y = locus_df[value]

            ax.plot(
                x,
                y,
                alpha=0.75,
                c=color_map[cluster],
                marker=defaults.line_plot_marker,
                markersize=defaults.line_plot_marker_size
            )

    ax.set_xlabel(defaults.sample_label, fontsize=defaults.axis_label_font_size)

    ax.set_ylabel(value.upper(), fontsize=defaults.axis_label_font_size)

    ax.set_xticks(sorted(df['sample_index'].unique()))

    ax.set_xticklabels(samples)

    utils.set_tick_label_font_sizes(ax, defaults.tick_label_font_size)

    ax.set_ylim(*defaults.cellular_prevalence_limits)

    box = ax.get_position()

    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    legend_handles = utils.get_legend_handles(color_map)

    legend = ax.legend(
        list(legend_handles.values()),
        list(legend_handles.keys()),
        bbox_to_anchor=(1, 0.5),
        fontsize=defaults.legend_font_size,
        loc='center left',
        title=defaults.cluster_label
    )

    legend.get_title().set_fontsize(defaults.legend_title_font_size)

    fig.set_size_inches(*utils.get_parallel_coordinates_figure_size(samples))

    utils.save_figure(fig, out_file)

#=======================================================================================================================
# Scatter
#=======================================================================================================================


def scatter_plot(
        config,
        trace,
        out_file,
        burnin=0,
        max_clusters=None,
        min_cluster_size=0,
        samples=None,
        thin=1,
        value='cellular_prevalence'):

    utils.setup_plot()

    df = post_process.loci.load_table(
        config,
        trace,
        burnin=burnin,
        max_clusters=max_clusters,
        min_cluster_size=min_cluster_size,
        thin=thin,
    )

    if len(samples) == 0:
        samples = sorted(df['sample'].unique())

    color_map = utils.get_clusters_color_map(df['cluster'])

    cluster_df = df[['mutation', 'cluster']].drop_duplicates().set_index('mutation')

    loci_color_map = cluster_df['cluster'].map(color_map).to_dict()

    mean_df = df.pivot(index='mutation', columns='sample', values=value)

    scatter.plot_all_pairs(
        loci_color_map,
        mean_df,
        out_file,
        samples,
        legend_color_map=color_map
    )

#=======================================================================================================================
# Similarity matrix
#=======================================================================================================================


def similarity_matrix_plot(
        trace,
        out_file,
        burnin=0,
        max_clusters=None,
        min_cluster_size=0,
        thin=1):

    sb.set_style('whitegrid')

    labels = post_process.cluster_pyclone_trace(trace, burnin, thin, max_clusters=max_clusters)

    labels = labels.set_index('mutation')

    labels = labels['cluster']

    color_map = utils.get_clusters_color_map(labels)

    cluster_sizes = labels.value_counts()

    used_clusters = cluster_sizes[cluster_sizes >= min_cluster_size].index

    labels = labels[labels.isin(used_clusters)]

    used_loci = labels.index

    labels_trace = trace.labels[used_loci]

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

    g.fig.savefig(out_file, bbox_inches='tight')

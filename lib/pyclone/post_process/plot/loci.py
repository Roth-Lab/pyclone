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

from pyclone.trace import load_cluster_labels_trace
from pyclone.post_process import cluster_pyclone_trace

import pyclone.paths as paths
import pyclone.post_process as post_process
import pyclone.trace as trace

import defaults
import utils


#=======================================================================================================================
# Density
#=======================================================================================================================
def density_plot(
    config_file,
    plot_file,
    axis_label_font_size=defaults.axis_label_font_size,
    burnin=0,
    min_cluster_size=0,
    samples=None,
    thin=1,
    tick_font_size=defaults.tick_font_size):
    
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
            fontsize=axis_label_font_size
        )
        
        utils.set_tick_label_font_sizes(ax, tick_font_size)

    fig.text(
        -0.01, 
        0.5, 
        defaults.cellular_prevalence_label,
        fontsize=axis_label_font_size,
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
# Scatter
#=======================================================================================================================
def scatter_plot(
    config_file,
    plot_file,
    axis_label_font_size=defaults.axis_label_font_size,
    burnin=0,
    min_cluster_size=0,
    samples=None,
    thin=1,
    tick_font_size=defaults.tick_font_size,
    value='cellular_prevalence'):
    
    utils.setup_plot()
    
    df = post_process.load_multi_sample_table(
        config_file,
        burnin,
        thin,
        min_cluster_size=min_cluster_size
    )
    
    if samples is None:
        samples = sorted(df['sample_id'].unique())
    
    else:
        df = df[df['sample_id'].isin(samples)]
    
    color_map = utils.get_clusters_color_map(df['cluster_id'])
    
    num_samples = len(samples)
    
    cluster_df = df[['mutation_id', 'cluster_id']].drop_duplicates()
    
    mean_df = df.pivot(index='mutation_id', columns='sample_id', values=value)
    
    size = max(2 * num_samples, 4)
    
    fig = pp.figure(figsize=(size, size))
    
    grid = gs.GridSpec(nrows=num_samples, ncols=num_samples)
    
    axes = []
    
    for i in range(num_samples):
        for j in range(i):
            ax = fig.add_subplot(grid[i, j])
    
            axes.append(ax)
    
            utils.setup_axes(ax)
            
            x_sample = samples[j]
            
            y_sample = samples[i]
            
            for cluster_id in color_map:
                loci = cluster_df.loc[cluster_df['cluster_id'] == cluster_id, 'mutation_id']
                
                x = mean_df.loc[loci][x_sample].values
        
                y = mean_df.loc[loci][y_sample].values
        
                ax.scatter(
                    x,
                    y,
                    s=20,
                    c=color_map[cluster_id],
                    label=cluster_id
                )
            
            ax.set_xlim(defaults.cellular_prevalence_limits)
            
            ax.set_ylim(defaults.cellular_prevalence_limits)
            
            utils.set_tick_label_font_sizes(ax, tick_font_size)
            
            if i != (num_samples - 1):
                ax.set_xticklabels([])
                
                ax.spines['bottom'].set_visible(False)
                
                ax.xaxis.set_ticks_position('none')
                
            else:
                ax.set_xlabel(x_sample, fontsize=axis_label_font_size)
                
            if j != 0:
                ax.set_yticklabels([])
            
                ax.spines['left'].set_visible(False)
                
                ax.yaxis.set_ticks_position('none')
                
            else:
                ax.set_ylabel(y_sample, fontsize=axis_label_font_size)

    legend = axes[0].legend(
        bbox_to_anchor=(1.1, 0.5),
        fontsize=6,
        loc='center left',
        title='Cluster'
    )
    
    legend.get_title().set_fontsize(8)
    
    utils.save_figure(fig, plot_file)

#=======================================================================================================================
# Similarity matrix
#=======================================================================================================================
def similarity_matrix_plot(config_file, plot_file, burnin=0, min_cluster_size=0, thin=1):
    
    sb.set_style('whitegrid')
    
    labels = cluster_pyclone_trace(config_file, burnin, thin)
    
    labels = labels.set_index('mutation_id')
    
    labels = labels['cluster_id']
    
    color_map = utils.get_clusters_color_map(labels)
    
    cluster_sizes = labels.value_counts()
    
    used_clusters = cluster_sizes[cluster_sizes >= min_cluster_size].index
    
    labels = labels[labels.isin(used_clusters)]

    used_loci = labels.index
          
    trace_file = paths.get_labels_trace_file(config_file)
    
    trace = load_cluster_labels_trace(trace_file, burnin, thin)

    dist_mat = pdist(trace.values.T, 'hamming')
    
    dist_mat = pd.DataFrame(squareform(dist_mat), index=trace.columns, columns=trace.columns)
    
    dist_mat = dist_mat.loc[used_loci, used_loci]
    
    sim_mat = 1 - dist_mat
   
    Z = average(dist_mat)
    
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
    
    utils.set_tick_label_font_sizes(ax, defaults.tick_font_size)
    
    utils.set_tick_label_rotations(ax)
            
    g.fig.savefig(plot_file, bbox_inches='tight')
    

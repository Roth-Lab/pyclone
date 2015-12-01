'''
Created on Nov 30, 2015

@author: Andrew Roth
'''
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D

import matplotlib.gridspec as gs
import matplotlib.pyplot as pp
import numpy as np
import pandas as pd
import seaborn as sb

from pyclone.post_process import load_cluster_posteriors_table, load_cluster_posteriors_summary_table

import defaults
import utils

def density_plot(
    config_file, 
    plot_file, 
    axis_label_font_size=defaults.axis_label_font_size,
    burnin=0, 
    mesh_size=101,
    min_cluster_size=0,
    samples=None,
    thin=1,
    tick_font_size=defaults.tick_font_size):
    
    df = load_cluster_posteriors_table(
        config_file, 
        burnin=burnin, 
        thin=thin, 
        mesh_size=mesh_size,
        min_size=min_cluster_size
    )

    sizes = df[['cluster_id', 'size']].drop_duplicates().set_index('cluster_id').to_dict()['size']
    
    if samples is None:
        samples = sorted(df['sample_id'].unique())
    
    num_samples = len(samples)
    
    clusters = df['cluster_id'].unique()
    
    postions = range(1, len(clusters) + 1)
    
    utils.setup_plot()
    
    width = 8
    
    height = 2 * num_samples + 1

    fig = pp.figure(figsize=(width, height))
    
    grid = gs.GridSpec(nrows=num_samples, ncols=1)
    
    colors = utils.get_clusters_color_map(pd.Series(clusters))
        
    for ax_index, sample_id in enumerate(samples):
        plot_df = df[df['sample_id'] == sample_id]
        
        plot_df = plot_df.drop(['sample_id', 'size'], axis=1).set_index('cluster_id')
        
        ax = fig.add_subplot(grid[ax_index])
                
        utils.setup_axes(ax)
    
        ax.annotate(
            sample_id,
            xy=(1.01, 0.5),
            xycoords='axes fraction', 
            fontsize=axis_label_font_size
        )
    
        for i, (cluster_id, log_pdf) in enumerate(plot_df.iterrows()):
            pos = postions[i]
    
            y = log_pdf.index.astype(float)
    
            x = np.exp(log_pdf)
            
            x = (x / x.max()) * 0.3
    
            ax.fill_betweenx(y, pos - x, pos + x, color=colors[cluster_id], where=(x>1e-6))
        
        ax.set_xticks(postions)
    
        if ax_index == (num_samples - 1):
            x_tick_labels = ['{0} (n={1})'.format(x, sizes[x]) for x in clusters]
            
            ax.set_xticklabels(
                x_tick_labels,
                rotation=90
            )
            
            ax.set_xlabel(defaults.cluster_label, fontsize=axis_label_font_size)
        
        else:
            ax.set_xticklabels([])
        
        utils.set_tick_label_font_sizes(ax, tick_font_size)
        
        ax.set_ylim(defaults.cellular_prevalence_limits)
    
    if num_samples == 1:
        ax.set_ylabel(
            defaults.cellular_prevalence_label,
            fontsize=axis_label_font_size
        )    
    
    else:
        fig.text(
            -0.01, 
            0.5, 
            defaults.cellular_prevalence_label,
            fontsize=axis_label_font_size,
            ha='center',
            rotation=90,
            va='center'
        )
        
    grid.tight_layout(fig)
    
    utils.save_figure(fig, plot_file)
    
def parallel_coordinates_plot(
    config_file,
    plot_file,
    axis_label_font_size=defaults.axis_label_font_size,
    burnin=0,
    mesh_size=101,
    min_cluster_size=0,
    samples=None,
    thin=1,
    tick_font_size=defaults.tick_font_size):
    
    utils.setup_plot()
    
    plot_df = load_cluster_posteriors_summary_table(
        config_file, 
        burnin=burnin,
        mesh_size=mesh_size,
        min_size=min_cluster_size,
        thin=thin, 
    )   
    
    if samples is None:
        samples = sorted(plot_df['sample_id'].unique())
    
    else:
        plot_df = plot_df[plot_df['sample_id'].isin(samples)]
    
    clusters = sorted(plot_df['cluster_id'].unique())
    
    plot_df['sample_index'] = plot_df['sample_id'].apply(lambda x: samples.index(x))
    
    plot_df = plot_df.sort_values(by='sample_index')
    
    grid = sb.FacetGrid(
        plot_df,
        hue='cluster_id',
        hue_order=clusters,
        palette='husl'
    )
    
    grid.map(
        pp.errorbar, 
        'sample_index', 
        'mean', 
        'std', 
        marker=defaults.line_plot_marker,
        markersize=defaults.line_plot_marker_size
    )
    
    ax = grid.ax
    
    utils.setup_axes(ax)
    
    fig = grid.fig
    
    # Legend
    box = ax.get_position()
    
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Cluster')
    
    # Axis formatting
    ax.set_xticks(sorted(plot_df['sample_index'].unique()))
    
    ax.set_xticklabels(samples)
    
    ax.set_xlabel(defaults.sample_label, fontsize=axis_label_font_size)
    
    ax.set_ylabel(defaults.cellular_prevalence_label, fontsize=axis_label_font_size)
    
    utils.set_tick_label_font_sizes(ax, tick_font_size)
    
    # Plot limits
    ax.set_xlim(
        plot_df['sample_index'].min() - 0.1,
        plot_df['sample_index'].max() + 0.1
    )
    
    ax.set_ylim(*defaults.cellular_prevalence_limits)
    
    # Resize and save figure    
    fig.set_size_inches(*utils.get_parallel_coordinates_figure_size(samples))
    
    utils.save_figure(fig, plot_file)
    
def scatter_plot(
    config_file,
    plot_file,
    axis_label_font_size=defaults.axis_label_font_size,
    burnin=0,
    mesh_size=101,
    min_cluster_size=0,
    samples=None,
    thin=1,
    tick_font_size=defaults.tick_font_size):
    
    utils.setup_plot()
    
    df = load_cluster_posteriors_summary_table(
        config_file, 
        burnin=burnin,
        mesh_size=mesh_size,
        min_size=min_cluster_size,
        thin=thin, 
    )   
    
    if samples is None:
        samples = sorted(df['sample_id'].unique())
    
    else:
        df = df[df['sample_id'].isin(samples)]
        
    color_map = utils.get_clusters_color_map(pd.Series(df['cluster_id']))
    
    num_samples = len(samples)
    
    grid = gs.GridSpec(nrows=num_samples, ncols=num_samples)
    
    size = max(2 * num_samples, 4)
    
    fig = pp.figure(figsize=(size, size))
    
    axes = []

    for i in range(num_samples):
        for j in range(i):
            ax = fig.add_subplot(grid[i, j])
            
            axes.append(ax)
            
            target_samples = [samples[j], samples[i]]
            
            plot_df = df[df['sample_id'].isin(target_samples)]
    
            _scatter_plot(ax, color_map, plot_df, target_samples[0], target_samples[1])
            
            utils.setup_axes(ax)
            
            ax.set_xlim(defaults.cellular_prevalence_limits)
            
            ax.set_ylim(defaults.cellular_prevalence_limits)

            utils.set_tick_label_font_sizes(ax, tick_font_size)
            
            if i != (num_samples - 1):
                ax.set_xticklabels([])
                
                ax.spines['bottom'].set_visible(False)
                
                ax.xaxis.set_ticks_position('none')
                
            else:
                ax.set_xlabel(target_samples[0], fontsize=axis_label_font_size)
                
            if j != 0:
                ax.set_yticklabels([])
            
                ax.spines['left'].set_visible(False)
                
                ax.yaxis.set_ticks_position('none')
                
            else:
                ax.set_ylabel(target_samples[1], fontsize=axis_label_font_size)
    
    legend_handles = []

    for cluster_id in color_map:
        legend_handles.append(
            Line2D(
                [0], 
                [0], 
                linestyle="none", 
                marker="o", 
                markersize=4, 
                markerfacecolor=color_map[cluster_id]
            )
        )
    
    legend = axes[0].legend(
        legend_handles, 
        [str(x) for x in color_map],
        bbox_to_anchor=(1.1, 0.5),
        fontsize=6,
        loc='center left',
        title='Cluster'
    )
        
    legend.get_title().set_fontsize(8)
 
    grid.tight_layout(fig)
    
    utils.save_figure(fig, plot_file)
    
def _scatter_plot(ax, color_map, df, x_sample, y_sample):
    coords = df.pivot(index='cluster_id', columns='sample_id', values='mean')

    error_bars = df.pivot(index='cluster_id', columns='sample_id', values='std')
    
    error_bars = error_bars.loc[coords.index]
    
    colors = [color_map[x] for x in coords.index]
    
    x = coords[x_sample].values
    
    y = coords[y_sample].values
    
    x_err = error_bars[x_sample].values
    
    y_err = error_bars[y_sample].values
   
    ax.scatter(x, y, alpha=0.8, c=colors, s=15)
    
    for i in range(len(x)):
        e = Ellipse(
            (x[i], y[i]),
            width=x_err[i],
            height=y_err[i],
        )
        
        e.set_facecolor(colors[i])
        
        e.set_alpha(0.2)
        
        ax.add_artist(e)
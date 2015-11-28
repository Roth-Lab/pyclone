'''
Created on 2013-06-06

@author: Andrew Roth
'''
from __future__ import division

import matplotlib.pyplot as pp
import seaborn as sb

from pyclone.post_process import load_multi_sample_table

from .utils import setup_axes, setup_plot

import pyclone.paths as paths

def plot_multi_sample_parallel_coordinates(config_file, plot_file, y_value, burnin=0, thin=1, samples=None, separate_lines=False):
    setup_plot()
    
    if samples is None:
        samples = paths.get_sample_ids(config_file)
    
    df = load_multi_sample_table(config_file, burnin, thin)
    
    if separate_lines:
        plot_df = df
        
    else:
        plot_df = df.groupby(['sample', 'cluster_id']).mean().reset_index()
    
    plot_df['sample_index'] = plot_df['sample'].apply(lambda x: samples.index(x))
    
    plot_df = plot_df.sort_values(by='sample_index')

    xticks = range(len(samples))
    
    if separate_lines:
        fig = pp.figure()
        
        ax = fig.add_subplot(1, 1, 1)
        
        sb.tsplot(
            ax=ax,
            
            data=df,
            condition='cluster_id',
            unit='mutation_id',
            time='sample_index', 
            value=y_value,
            
            err_style='unit_traces',
            marker="o",
            markersize=4
        )
        
    else:
        grid = sb.FacetGrid(
            plot_df,
            hue='cluster_id',
            hue_order=sorted(plot_df['cluster_id'].unique()),
            palette='husl'
        )
        
        ax = grid.ax
        
        fig = grid.fig 
        
        if y_value == 'cellular_prevalence':
            grid.map(
                pp.errorbar, 
                'sample_index', 
                'cellular_prevalence', 
                'cellular_prevalence_std', 
                marker="o",
                markersize=4
            )
        
        elif y_value == 'variant_allele_frequency':
            grid.map(
                pp.plot, 
                'sample_index', 
                'variant_allele_frequency', 
                marker="o",
                markersize=4
            )

    setup_axes(ax)
        
    ax.set_xlim(min(xticks) - 0.1, max(xticks) + 0.1)
    
    ax.set_ylim(0, 1)
    
    ax.set_xlabel('Sample')
    
    if y_value == 'cellular_prevalence':
        ax.set_ylabel('Cellular prevalence')
    
    elif y_value == 'variant_allele_frequency':
        ax.set_ylabel('VAF')
    
    ax.set_xticks(xticks)
    
    ax.set_xticklabels(samples, size=8, rotation=90)
    
    # Shrink current axis by 20%
    box = ax.get_position()
    
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Cluster')

    height = 4
    
    width = 1 * len(samples) + 2
        
    fig.set_figheight(height)
    
    fig.set_figwidth(width)
    
    fig.savefig(plot_file, bbox_inches='tight')

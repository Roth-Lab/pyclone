'''
Created on Dec 2, 2015

@author: Andrew Roth
'''
from __future__ import division

from math import ceil
from matplotlib.patches import Ellipse

import matplotlib.gridspec as gs
import matplotlib.pyplot as pp

import defaults
import utils

def plot_all_pairs(color_map, mean_df, plot_file, samples, error_df=None, legend_color_map=None):
    
    num_samples = len(samples)
    
    grid = gs.GridSpec(nrows=num_samples, ncols=num_samples)
    
    size = max(2 * num_samples, 4)
    
    fig = pp.figure(figsize=(size, size))
    
    axes = []

    for i in range(num_samples):
        for j in range(i):
            ax = fig.add_subplot(grid[i, j])
            
            axes.append(ax)
            
            x_sample = samples[j]
            
            y_sample = samples[i]
    
            _plot(
                ax, 
                color_map, 
                mean_df, 
                x_sample, 
                y_sample,
                error_df=error_df
            )
            
            utils.setup_axes(ax)

            utils.set_tick_label_font_sizes(ax, defaults.tick_label_font_size)
            
            if i != (num_samples - 1):
                ax.set_xticklabels([])
                
                ax.spines['bottom'].set_visible(False)
                
                ax.xaxis.set_ticks_position('none')
                
            else:
                ax.set_xlabel(x_sample, fontsize=defaults.axis_label_font_size)
                
            if j != 0:
                ax.set_yticklabels([])
            
                ax.spines['left'].set_visible(False)
                
                ax.yaxis.set_ticks_position('none')
                
            else:
                ax.set_ylabel(y_sample, fontsize=defaults.axis_label_font_size)
    
    if legend_color_map is None:
        legend_color_map = color_map
    
    add_legend(axes[0], legend_color_map)
   
    grid.tight_layout(fig)
    
    utils.save_figure(fig, plot_file)

def _plot(ax, color_map, df, x_sample, y_sample, error_df=None):
    
    colors = [color_map[x] for x in df.index]
    
    x = df[x_sample].values
    
    y = df[y_sample].values
        
    if error_df is not None:
        error_df = error_df.loc[df.index]
    
        x_err = error_df[x_sample].values
    
        y_err = error_df[y_sample].values
   
    ax.scatter(x, y, alpha=0.8, c=colors, s=15)
    
    for i in range(len(x)):
        if error_df is not None:
            e = Ellipse(
                (x[i], y[i]),
                width=x_err[i],
                height=y_err[i],
            )
            
            e.set_facecolor(colors[i])
            
            e.set_alpha(0.2)
            
            ax.add_artist(e)
    
    ax.set_xlim(*defaults.cellular_prevalence_limits)
    
    ax.set_ylim(*defaults.cellular_prevalence_limits)
    
def add_legend(ax, color_map):
    
    legend_handles = utils.get_legend_handles(
        color_map, 
        marker=defaults.line_plot_marker, 
        markersize=defaults.line_plot_marker_size
    )
    
    num_cols = int(ceil(len(legend_handles) / 8))
    
    legend = ax.legend(
        legend_handles.values(), 
        legend_handles.keys(),
        bbox_to_anchor=(1.1, 0.5),
        fontsize=defaults.legend_font_size,
        loc='center left',
        ncol=num_cols,
        title='Cluster'
    )
        
    legend.get_title().set_fontsize(defaults.legend_title_font_size)
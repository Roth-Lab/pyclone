'''
Created on 2012-08-20

@author: Andrew Roth
'''
from collections import OrderedDict
from matplotlib.lines import Line2D

import seaborn as sb

def get_legend_handles(color_map, marker='o', markersize=4):
    legend_handles = OrderedDict()

    for key in color_map:
        legend_handles[str(key)] = Line2D(
            [0],
            [0],
            linestyle='none',
            marker=marker,
            markersize=markersize,
            markerfacecolor=color_map[key]
        )
    
    return legend_handles

def get_clusters_color_map(labels, palette='husl'):
    clusters = sorted(labels.unique())
    
    num_clusters = len(clusters)
        
    color_map = OrderedDict(zip(clusters, sb.color_palette(palette, num_clusters)))
    
    return color_map

def get_parallel_coordinates_figure_size(samples):
    width = 1 * len(samples) + 2
    
    height = 4
    
    return width, height

def save_figure(fig, file_name):
    fig.savefig(file_name, bbox_inches='tight')
    
def set_axis_label_font_size(ax, size):
    ax.set_xlabel(ax.get_xlabel(), fontsize=size)
    
    ax.set_ylabel(ax.get_ylabel(), fontsize=size)

def set_tick_label_font_sizes(ax, size):
    for t in ax.get_xticklabels():
        t.set_size(size)
    
    for t in ax.get_yticklabels():
        t.set_size(size)

def set_tick_label_rotations(ax):
    for t in ax.get_xticklabels():
        t.set_rotation(90)
        
    for t in ax.get_yticklabels():
        t.set_rotation(0)

def setup_plot():    
    sb.set_style('ticks', {'font.sans-serif':['Helvetica']})

def setup_axes(ax):
    ax.spines['left'].set_position(('outward', 10))
    
    ax.spines['bottom'].set_position(('outward', 10))
    
    ax.spines['top'].set_visible(False)
    
    ax.spines['right'].set_visible(False)
    
    ax.xaxis.tick_bottom()
    
    ax.yaxis.tick_left()
    
    ax.xaxis.grid(True, which="major", linestyle=':')
    
    ax.yaxis.grid(True, which="major", linestyle=':')

    return ax

'''
Created on 2012-08-20

@author: Andrew Roth
'''
import seaborn as sb

def get_clusters_color_map(labels, palette='husl'):
    clusters = sorted(labels.unique())
    
    num_clusters = len(clusters)
        
    color_map = dict(zip(clusters, sb.color_palette(palette, num_clusters)))
    
    return color_map

def get_parallel_coordinates_figure_size(samples):
    width = 1 * len(samples) + 2
    
    height = 4
    
    return width, height

def save_figure(fig, file_name):
    fig.savefig(file_name, bbox_inches='tight')

def set_tick_font_sizes(ax, size):
    for t in ax.get_xticklabels():
        t.set_size(size)
    
    for t in ax.get_yticklabels():
        t.set_size(size)

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

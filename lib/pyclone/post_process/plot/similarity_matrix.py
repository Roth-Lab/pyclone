'''
Created on Nov 27, 2015

@author: Andrew Roth
'''
from scipy.cluster.hierarchy import average
from scipy.spatial.distance import pdist, squareform

import pandas as pd
import seaborn as sb

from pyclone.post_process import cluster_pyclone_trace
from pyclone.trace import load_cluster_labels_trace

from .utils import get_clusters_color_map

import pyclone.paths as paths

def plot_similarity_matrix(config_file, plot_file, burnin, thin):
    sb.set_style("whitegrid")
    
    trace_file = paths.get_labels_trace_file(config_file)
    
    trace = load_cluster_labels_trace(trace_file, burnin, thin)

    dist_mat = pdist(trace.values.T, 'hamming')
    
    dist_mat = squareform(dist_mat)
    
    sim_mat = 1 - dist_mat
    
    sim_mat = pd.DataFrame(sim_mat, index=trace.columns, columns=trace.columns) 
    
    Z = average(dist_mat)
    
    N = sim_mat.shape[0]
    
    labels = cluster_pyclone_trace(config_file, burnin, thin)
    
    labels = labels.set_index('mutation_id')
    
    labels = labels.loc[sim_mat.index]['cluster_id']
    
    cluster_color_map = get_clusters_color_map(labels)
    
    # Convert the palette to vectors that will be drawn on the side of the matrix
    cluster_colors = labels.map(cluster_color_map)
    
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

    for t in ax.get_xticklabels():
        t.set_rotation(90)
        
        t.set_size(6)
        
    for t in ax.get_yticklabels():
        t.set_rotation(0)
        
        t.set_size(6)
            
    g.fig.savefig(plot_file, bbox_inches='tight')
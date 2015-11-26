from scipy.cluster.hierarchy import average
from scipy.spatial.distance import pdist, squareform

import pandas as pd
import seaborn as sb

from pyclone.post_process.utils import load_cellular_frequencies_trace, load_cluster_labels_trace
from pyclone.post_process.cluster import cluster_pyclone_trace

def plot_similarity_matrix(trace_file, plot_file, burnin, thin):
    trace = load_cluster_labels_trace(trace_file, burnin, thin)

    dist_mat = pdist(trace.values.T, 'hamming')
    
    dist_mat = squareform(dist_mat)
    
    sim_mat = 1 - dist_mat
    
    sim_mat = pd.DataFrame(sim_mat, index=trace.columns, columns=trace.columns) 
    
    Z = average(dist_mat)
    
    N = sim_mat.shape[0]
    
    labels = cluster_pyclone_trace(trace_file, burnin, thin)
    
    labels = labels.set_index('mutation_id')
    
    labels = labels.loc[sim_mat.index]['cluster_id']
    
    used_clusters = labels.unique()
    
    # Create a custom palette to identify the networks
    cluster_palette = sb.color_palette('husl', len(used_clusters))
  
    cluster_color_map = dict(zip(used_clusters, cluster_palette))
    
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

def plot_cellular_frequencies(trace_file, plot_file, burnin, thin):
    pass
#     trace = load_cellular_frequencies_trace(trace_file, burnin, thin)
#         
#     plotter = CellularFrequencyPlot(trace)
# 
#     plotter.plot()
#     
#     plotter.save(plot_file)

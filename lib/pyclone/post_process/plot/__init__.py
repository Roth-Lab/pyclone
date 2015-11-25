from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import average

import pandas as pd
import seaborn as sb

def plot_similarity_matrix(trace, out_file):
    d = squareform(pdist(trace.T, metric='hamming'))
    
    s = 1 - d
    
    z = average(d)
    
    s = pd.DataFrame(s, index=trace.columns, columns=trace.columns)
    
    dim = s.shape[0] * 0.16

    g = sb.clustermap(s,
                      figsize=(dim, dim),
                      col_linkage=z, 
                      row_linkage=z)
    
    for t in g.ax_heatmap.get_xticklabels():
        t.set_fontsize(8)
        
    for t in g.ax_heatmap.get_yticklabels():
        t.set_fontsize(8)
        
    g.savefig(out_file, bbox_inches='tight')

# def plot_cellular_frequencies(trace_file, plot_file, burnin, thin):
#     trace = load_cellular_frequencies_trace(trace_file, burnin, thin)
#         
#     plotter = CellularFrequencyPlot(trace, cmap=colors)
# 
#     plotter.plot()
#     
#     plotter.save(plot_file)

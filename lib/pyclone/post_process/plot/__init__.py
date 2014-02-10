from pyclone.post_process.utils import load_cellular_frequencies_trace, load_cluster_labels_trace

from cellular_frequencies import CellularFrequencyPlot
from densities import PosteriorDensity
from similarity_matrix import SimilarityMatrixPlot

import brewer2mpl

bmap = brewer2mpl.get_map('Set1', 'qualitative', 9)

colors = bmap.mpl_colormap

def plot_similarity_matrix(trace_file, plot_file, burnin, thin):
    trace = load_cluster_labels_trace(trace_file, burnin, thin)
    
    if len(trace) < 2:
        return
    
    plotter = SimilarityMatrixPlot(trace)
    
    plotter.plot()
    
    plotter.save(plot_file)

def plot_cellular_frequencies(trace_file, plot_file, burnin, thin):
    trace = load_cellular_frequencies_trace(trace_file, burnin, thin)
        
    plotter = CellularFrequencyPlot(trace, cmap=colors)

    plotter.plot()
    
    plotter.save(plot_file)

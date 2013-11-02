from pyclone.post_process.utils import load_trace

from cellular_frequencies import CellularFrequencyPlot
from densities import PosteriorDensity
from similarity_matrix import SimilarityMatrixPlot

def plot_similarity_matrix(trace_file, plot_file, burnin, thin):
    trace = load_trace(trace_file, burnin, thin)
    
    if len(trace) < 2:
        return
    
    plotter = SimilarityMatrixPlot(trace)
    
    plotter.plot()
    
    plotter.save(plot_file)

def plot_cellular_frequencies(trace_file, plot_file, burnin, thin):
    trace = load_trace(trace_file, burnin, thin)
        
    plotter = CellularFrequencyPlot(trace)

    plotter.plot()
    
    plotter.save(plot_file)

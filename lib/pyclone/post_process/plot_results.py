import bz2
import brewer2mpl
import csv
import numpy as np
import os

from collections import defaultdict
from pyclone.plot import CellularFrequencyPlot, SimilarityMatrixPlot
from scipy.spatial.distance import pdist, squareform

csv.field_size_limit(10000000)

bmap = brewer2mpl.get_map('Set1', 'qualitative', 9)

colors = bmap.mpl_colormap

def main(trace_file, plot_file, plot_type, burnin, thin):    
    if plot_type == 'similarity_matrix':
        plot_similarity_matrix(trace_file, plot_file, burnin, thin)

    elif plot_type == 'cellular_frequencies':
        plot_cellular_frequencies(trace_file, plot_file, burnin, thin)

def plot_similarity_matrix(trace_file, plot_file, burnin, thin):
    trace = load_trace(trace_file, burnin, thin, int)
    
    if len(trace) < 2:
        return
    
    plotter = SimilarityMatrixPlot(trace)
    
    plotter.plot()
    
    plotter.save(plot_file)

def plot_cellular_frequencies(trace_file, plot_file, burnin, thin):
    trace = load_trace(trace_file, burnin, thin, float)
        
    plotter = CellularFrequencyPlot(trace, cmap=colors)

    plotter.plot()
    
    plotter.save(plot_file)

def load_trace(trace_file, burnin, thin, cast_func):
    '''
        Args:
            trace_file : (str) Path to file to load.
            burnin : (int) Number of samples from the begining of MCMC chain to discard.
            thin : (int) Number of samples to skip when building trace.
            cast_func : (function) A function to cast data from string to appropriate type i.e. int, float
    '''        
    trace = defaultdict(list)
    
    fh = bz2.BZ2File(trace_file)
    
    reader = csv.DictReader(fh, delimiter='\t')
    
    for i, row in enumerate(reader):
        if i < burnin:
            continue
        
        if i % thin == 0:
            for mutation in row:
                trace[mutation].append(cast_func(row[mutation]))
        
    fh.close()
    
    return trace

if __name__ == "__main__":
    import sys
        
    trace_file = sys.argv[1]
    
    plot_file = sys.argv[2]
    
    plot_type = sys.argv[3]
    
    burnin = int(sys.argv[4])
    
    thin = int(sys.argv[5])
            
    main(trace_file, plot_file, plot_type, burnin, thin)

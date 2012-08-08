from pyclone.plot import CellularFrequencyPlot, SimilarityMatrixPlot
from pyclone.post_process import DpSamplerPostProcessor
from pyclone.results import AnalysisDB

import matplotlib.pyplot as plot
import os

def main(pyclone_file, out_folder):
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    
    results = AnalysisDB(pyclone_file)
    
    post_processor = DpSamplerPostProcessor(results)
    
    sm = SimilarityMatrixPlot()
    
    genes = post_processor.genes
    
    sim_mat = post_processor.get_similarity_posteriors(1000, 10)
    
    sm.plot(genes, sim_mat)
    
    sm_file = os.path.join(out_folder, 'similarity_matrix.pdf')
    
    sm.save(sm_file)
    
    clusters = sm.get_clusters(genes, sim_mat)
        
    cellular_frequencies = post_processor.cellular_frequencies
    
    cp = CellularFrequencyPlot(cellular_frequencies, clusters)
    
    cp.plot()
    
    cp_file = os.path.join(out_folder, 'cellular_frequencies.pdf')
    
    cp.save(cp_file)
    
    results.close()

if __name__ == "__main__":
    import sys
    
    pyclone_file = sys.argv[1]
    out_folder = sys.argv[2]
    
    main(pyclone_file, out_folder)

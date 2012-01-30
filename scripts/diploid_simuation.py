'''
Simulation of single clonal population with diploid genomes with a mixture of heterozygous and homozygous mutations.
'''
from __future__ import division

from pyclone.simulation.data_simulator import ClonalDataSimulator

import matplotlib.pyplot as plot

simulator = ClonalDataSimulator(1, 2, 2, 1000)

allele_freqs = []

for _ in range(1000):
    a, d, cn, genotype, clone_freq = simulator.draw_data_point()
    
    allele_freqs.append(a / d) 

fig = plot.figure()
ax = fig.add_subplot(1, 1, 1)
ax.hist(allele_freqs, bins=100)
ax.axvline(clone_freq, c='y')
ax.set_xlim(0, 1)
plot.show()

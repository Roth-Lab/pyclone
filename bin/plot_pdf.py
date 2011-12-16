from __future__ import division

from matplotlib import rc

import csv
import glob
import os
import matplotlib.pyplot as plot

rc('text', usetex=True)
rc('font', family='serif')

data_dir = "../plots"

files = glob.glob(os.path.join(data_dir, "*.density.tsv"))

all_fig = plot.figure()
all_ax = plot.subplot(1, 1, 1)

all_ax.set_title(r'All Genes')
all_ax.set_xlabel(r'$\phi^i$')
all_ax.set_ylabel(r'$P(\phi^i)$')

for file_name in files:
    base_name = os.path.basename(file_name)
    
    parts = base_name.split('.')
    
    gene = parts[0]
    
    reader = csv.reader(open(file_name), delimiter='\t')
    
    x = reader.next()
    y = reader.next()

    plot_file = os.path.join(data_dir, '{0}.pdf'.format(gene))

    plot.figure()
    
    plot.plot(x, y)
    
    plot.title(r'{0}'.format(gene))
    plot.xlabel(r'$\phi^i$')
    plot.ylabel(r'$P(\phi^i)$')

    plot.savefig(plot_file)
    
    plot.close()
    
    all_ax.plot(x, y)

plot_file = os.path.join(data_dir, 'all.pdf')    
all_fig.savefig(plot_file)
    


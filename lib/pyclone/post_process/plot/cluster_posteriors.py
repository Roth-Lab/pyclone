'''
Created on Nov 30, 2015

@author: Andrew Roth
'''
import matplotlib.gridspec as gs
import matplotlib.pyplot as pp
import numpy as np
import pandas as pd

from pyclone.post_process import load_cluster_posteriors_table

from .utils import get_clusters_color_map, setup_axes, setup_plot

def plot_cluster_posteriors(
    config_file, 
    plot_file, 
    burnin, 
    thin,
    label_font_size=12,
    mesh_size=101,
    samples=None,
    sample_label_font_size=12,
    tick_font_size=8):
    
    df = load_cluster_posteriors_table(
        config_file, 
        burnin, 
        thin, 
        mesh_size=mesh_size
    )

    sizes = df[['cluster_id', 'size']].drop_duplicates().set_index('cluster_id').to_dict()['size']
    
    if samples is None:
        samples = sorted(df['sample_id'].unique())
    
    num_samples = len(samples)
    
    clusters = df['cluster_id'].unique()
    
    postions = range(1, len(clusters) + 1)
    
    setup_plot()
    
    width = 8
    
    height = 2 * num_samples + 1

    fig = pp.figure(figsize=(width, height))
    
    grid = gs.GridSpec(nrows=num_samples, ncols=1)
    
    colors = get_clusters_color_map(pd.Series(clusters))
        
    for ax_index, sample_id in enumerate(samples):
        plot_df = df[df['sample_id'] == sample_id]
        
        plot_df = plot_df.drop(['sample_id', 'size'], axis=1).set_index('cluster_id')
        
        ax = fig.add_subplot(grid[ax_index])
                
        setup_axes(ax)
    
        ax.annotate(
            sample_id,
            xy=(1.01, 0.5),
            xycoords='axes fraction', 
            fontsize=sample_label_font_size
        )
    
        for i, (cluster_id, log_pdf) in enumerate(plot_df.iterrows()):
            pos = postions[i]
    
            y = log_pdf.index.astype(float)
    
            x = np.exp(log_pdf)
            
            x = (x / x.max()) * 0.3
    
            ax.fill_betweenx(y, pos - x, pos + x, color=colors[cluster_id], where=(x>1e-6))
        
        ax.set_xticks(postions)
    
        if ax_index == (num_samples - 1):
            x_tick_labels = ['{0} (n={1})'.format(x, sizes[x]) for x in clusters]
            
            ax.set_xticklabels(
                x_tick_labels,
                fontsize=tick_font_size,
                rotation=90
            )
            
            ax.set_xlabel('Cluster', fontsize=label_font_size)
        
        else:
            ax.set_xticklabels([])
        
        for t in ax.get_yticklabels():
            t.set_size(tick_font_size)
        
        ax.set_ylim(-0.01, 1.01)
    
    if num_samples == 1:
        ax.set_ylabel(
            'Cellular prevalence',
            fontsize=label_font_size
        )    
    
    else:
        fig.text(
            -0.01, 
            0.5, 
            'Cellular prevalence',
            fontsize=label_font_size,
            ha='center',
            rotation=90,
            va='center'
        )
        
    grid.tight_layout(fig)
    
    fig.savefig(plot_file, bbox_inches='tight')
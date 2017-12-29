import click

import pyclone.run as run


#=======================================================================================================================
# Build tables
#=======================================================================================================================
@click.command(
    context_settings={'max_content_width': 120},
    name='table'
)
@click.option(
    '-t', '--trace-file',
    required=True,
    type=click.Path(exists=True, resolve_path=True),
    help='''Path to trace file will be written in HDF5 format.'''
)
@click.option(
    '-o', '--out-file',
    required=True,
    type=click.Path(resolve_path=True),
    help='''Path where table will be written in tab delimited format.'''
)
@click.option(
    '-f', '--table-format',
    required=True,
    type=click.Choice(['cluster', 'loci', 'old']),
    help='''Build a table of results. Choices are: `cluster` for cluster specific information; `loci` for loci specific
    information; `old` matches the 0.12.x PyClone output.'''
)
@click.option(
    '-m', '--max-clusters',
    default=100,
    type=int,
    help='''Maximum number of clusters to consider for post-processing.Note this does not affect the DP sampling only
    the final post-processing steps to get hard cluster assignments. Default is 100.'''
)
@click.option(
    '--burnin',
    default=0,
    type=int,
    help='''Number of samples to discard as burning for the MCMC chain. Default is 0.'''
)
@click.option(
    '--thin',
    default=1,
    type=int,
    help='''Number of samples to thin MCMC trace. For example if thin=10 every tenth sample after burning will be used
    for inference. Default is 1.'''
)
@click.option(
    '--grid-size',
    default=101,
    type=int,
    help='''Number of points to use for approximating the cluster posteriors. Most users should not need to change this.
    Default is 101.'''
)
def build_table(**kwargs):
    """ Build a summary table of a PyClone analysis.
    """
    run.build_table(**kwargs)


#=======================================================================================================================
# Plotting
#=======================================================================================================================
@click.command(
    context_settings={'max_content_width': 120},
    name='clusters'
)
@click.option(
    '-t', '--trace-file',
    required=True,
    type=click.Path(exists=True, resolve_path=True),
    help='''Path to trace file will be written in HDF5 format.'''
)
@click.option(
    '-o', '--out-file',
    required=True,
    type=click.Path(resolve_path=True),
    help='''Path to file where plot will be saved. Format can be controlled by changing file extension.'''
)
@click.option(
    '-f', '--plot-format',
    required=True,
    type=click.Choice(['density', 'line', 'scatter']),
    help='''Determines which style of plot will be done. Choices are: `denisty` plots posterior cluster ccf density,
    `line` parallel coordinate plots of mean ccf with error bars representing std and `scatter` plots grid of pairwise
    scatters of mean ccf with size proportional to std.'''
)
@click.option(
    '-m', '--max-clusters',
    default=100,
    type=int,
    help='''Maximum number of clusters to consider for post-processing.Note this does not affect the DP sampling only
    the final post-processing steps to get hard cluster assignments. Default is 100.'''
)
@click.option(
    '-s', '--samples',
    multiple=True,
    help='''Sample to plot. Can be specified multiple times for multiple samples. The order samples are specified
    controls the plotting order.'''
)
@click.option(
    '--burnin',
    default=0,
    type=int,
    help='''Number of samples to discard as burning for the MCMC chain. Default is 0.'''
)
@click.option(
    '--min-cluster-size',
    default=0,
    type=int,
    help='''Clusters with fewer mutations than this value will not be plotted. Default is to plot all clusters.'''
)
@click.option(
    '--thin',
    default=1,
    type=int,
    help='''Number of samples to thin MCMC trace. For example if thin=10 every tenth sample after burning will be used
    for inference. Default is 1.'''
)
@click.option(
    '--grid-size',
    default=101,
    type=int,
    help='''Number of points to use for approximating the cluster posteriors. Most users should not need to change this.
    Default is 101.'''
)
def plot_clusters(**kwargs):
    """ Plot results by cluster.
    """
    run.plot_clusters(**kwargs)


@click.command(
    context_settings={'max_content_width': 120},
    name='loci'
)
@click.option(
    '-t', '--trace-file',
    required=True,
    type=click.Path(exists=True, resolve_path=True),
    help='''Path to trace file will be written in HDF5 format.'''
)
@click.option(
    '-o', '--out-file',
    required=True,
    type=click.Path(resolve_path=True),
    help='''Path to file where plot will be saved. Format can be controlled by changing file extension.'''
)
@click.option(
    '-f', '--plot-format',
    required=True,
    type=click.Choice([
        'ccf-density', 'ccf-line', 'ccf-scatter', 'similarity-matrix', 'vaf-line', 'vaf-scatter'
    ]),
    help='''Determines which style of plot will be done.'''
)
@click.option(
    '-m', '--max-clusters',
    default=100,
    type=int,
    help='''Maximum number of clusters to consider for post-processing.Note this does not affect the DP sampling only
    the final post-processing steps to get hard cluster assignments. Default is 100.'''
)
@click.option(
    '-s', '--samples',
    multiple=True,
    help='''Sample to plot. Can be specified multiple times for multiple samples. The order samples are specified
    controls the plotting order.'''
)
@click.option(
    '--burnin',
    default=0,
    type=int,
    help='''Number of samples to discard as burning for the MCMC chain. Default is 0.'''
)
@click.option(
    '--min-cluster-size',
    default=0,
    type=int,
    help='''Clusters with fewer mutations than this value will not be plotted. Default is to plot all clusters.'''
)
@click.option(
    '--thin',
    default=1,
    type=int,
    help='''Number of samples to thin MCMC trace. For example if thin=10 every tenth sample after burning will be used
    for inference. Default is 1.'''
)
def plot_loci(**kwargs):
    """ Plot results by loci.
    """
    run.plot_loci(**kwargs)


#=======================================================================================================================
# Analysis
#=======================================================================================================================
@click.command(
    context_settings={'max_content_width': 120},
    name='new'
)
@click.option(
    '-i', '--in-file',
    required=True,
    type=click.Path(exists=True, resolve_path=True),
    help='''Path to TSV format file with copy number and allele count information for all samples. See the examples
    directory in the GitHub repository for format.'''
)
@click.option(
    '-t', '--trace-file',
    required=True,
    type=click.Path(resolve_path=True),
    help='''Path to where trace file will be written in HDF5 format.'''
)
@click.option(
    '-d', '--density',
    default='beta-binomial',
    type=click.Choice(['binomial', 'beta-binomial']),
    help='''Allele count density in the PyClone model. Use beta-binomial for most cases. Default beta-binomial.'''
)
@click.option(
    '-n', '--num-iters',
    default=int(1e4),
    type=int,
    help='''Number of iterations of the MCMC sampler to perform. Default is 10,000.'''
)
@click.option(
    '--concentration-value',
    default=1.0,
    type=float,
    help='''The (initial) concentration of the Dirichlet process. Higher values will encourage more clusters, lower
    values have the opposite effect. Default is 1.0.'''
)
@click.option(
    '--config-file',
    type=click.Path(exists=True, resolve_path=True),
    help='''Path to configuration file. Entries in this file will override the values set on the command line. See the
    examples directory in the GitHub repository for format.'''
)
@click.option(
    '--grid-size',
    default=None,
    type=int,
    help='''Grid size for discrete approximation. This will numerically marginalise the cancer cell fraction. Higher
    values lead to more accurate approximations at the expense of run time. By default this off and cancer cell
    fractions are sampled exactly. '''
)
@click.option(
    '--no-concentration-update',
    is_flag=True,
    help='''Set this to disable updating the concentration parameter of the Dirichlet process. Has not effect when the
    Binomial is used.'''
)
@click.option(
    '--no-precision-update',
    is_flag=True,
    help='''Set this to disable updating the precision parameter of the Beta-Binomial density. Has not effect when the
    Binomial is used.'''
)
@click.option(
    '--precision-value',
    default=400,
    type=float,
    help='''The (initial) precision parameter of the Beta-Binomial density. The higher the value the more similar the
    Beta-Binomial is to a Binomial. Default is 400.'''
)
@click.option(
    '--seed',
    default=None,
    type=int,
    help='''Set random seed so results can be reproduced. By default a random seed is chosen.'''
)
def new_analysis(**kwargs):
    """ Run a new PyClone analysis.
    """
    run.run_analysis(**kwargs)


@click.command(
    context_settings={'max_content_width': 120},
    name='resume'
)
@click.option(
    '-t', '--trace-file',
    required=True,
    type=click.Path(resolve_path=True),
    help='''Path to trace file will be written in HDF5 format.'''
)
@click.option(
    '-n', '--num-iters',
    default=int(1e4),
    type=int,
    help='''Number of iterations of the MCMC sampler to perform. Default is 10,000.'''
)
def resume_analysis(**kwargs):
    """ Resume a previous PyClone analysis from last MCMC iteration.
    """
    run.resume_analysis(**kwargs)


#=======================================================================================================================
# Setup main interface
#=======================================================================================================================
@click.group()
def pyclone():
    pass


@click.group()
def analysis():
    pass


@click.group(name='post-process')
def post_process():
    """ Post process results of a PyClone analysis.
    """
    pass


@click.group()
def plot():
    """ Plot the results of a PyClone analysis.
    """
    pass


pyclone.add_command(analysis)
pyclone.add_command(post_process)

analysis.add_command(new_analysis)
analysis.add_command(resume_analysis)

post_process.add_command(build_table)
post_process.add_command(plot)

plot.add_command(plot_clusters)
plot.add_command(plot_loci)

import click

import pyclone.run as run


#=======================================================================================================================
# Build tables
#=======================================================================================================================
@click.command(
    context_settings={'max_content_width': 120},
    name='build-table',
    help='Build a summary table of a PyClone analysis.'
)
@click.option(
    '-c', '--config-file',
    required=True,
    type=click.Path(exists=True, resolve_path=True),
    help='''Path to configuration file used for analysis. Use pyclone setup-analysis to build this file.'''
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
    run.build_table(**kwargs)


#=======================================================================================================================
# Plot clusters
#=======================================================================================================================
@click.command(
    context_settings={'max_content_width': 120},
    name='plot-clusters',
    help='Plot results by cluster.'
)
@click.option(
    '-c', '--config-file',
    required=True,
    type=click.Path(exists=True, resolve_path=True),
    help='''Path to configuration file used for analysis. Use pyclone setup-analysis to build this file.'''
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
    run.plot_clusters(**kwargs)


#=======================================================================================================================
# Plot loci
#=======================================================================================================================
@click.command(
    context_settings={'max_content_width': 120},
    name='plot-loci',
    help='Plot results by loci.'
)
@click.option(
    '-c', '--config-file',
    required=True,
    type=click.Path(exists=True, resolve_path=True),
    help='''Path to configuration file used for analysis. Use pyclone setup-analysis to build this file.'''
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
        'ccf_density', 'ccf_line', 'ccf_scatter', 'similarity_matrix', 'vaf_line', 'vaf_scatter'
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
    run.plot_loci(**kwargs)

#=======================================================================================================================
# Resume analysis
#=======================================================================================================================


@click.command(
    context_settings={'max_content_width': 120},
    name='resume-analysis',
    help='Resume a previous PyClone analysis.'
)
@click.option(
    '-c', '--config-file',
    required=True,
    type=click.Path(exists=True, resolve_path=True),
    help='''Path to configuration file used for analysis. Use pyclone setup-analysis to build this file.'''
)
@click.option(
    '-n', '--num-iters',
    default=int(1e4),
    type=int,
    help='''Number of iterations of the MCMC sampler to perform. Default is 10,000.'''
)
@click.option(
    '-t', '--trace-file',
    required=True,
    type=click.Path(resolve_path=True),
    help='''Path to trace file will be written in HDF5 format.'''
)
def resume_analysis(**kwargs):
    run.resume_analysis(**kwargs)


#=======================================================================================================================
# Run analysis
#=======================================================================================================================
@click.command(
    context_settings={'max_content_width': 120},
    name='run-analysis',
    help='Run an PyClone analysis.'
)
@click.option(
    '-c', '--config-file',
    required=True,
    type=click.Path(exists=True, resolve_path=True),
    help='''Path to configuration file used for analysis. Use pyclone setup-analysis to build this file.'''
)
@click.option(
    '-n', '--num-iters',
    default=int(1e4),
    type=int,
    help='''Number of iterations of the MCMC sampler to perform. Default is 10,000.'''
)
@click.option(
    '-t', '--trace-file',
    required=True,
    type=click.Path(resolve_path=True),
    help='''Path to where trace file will be written in HDF5 format.'''
)
@click.option(
    '--seed',
    default=None,
    type=int,
    help='''Set random seed so results can be reproduced. By default a random seed is chosen.'''
)
def run_analysis(**kwargs):
    run.run_analysis(**kwargs)


#=======================================================================================================================
# Setup analysis
#=======================================================================================================================
@click.command(
    context_settings={'max_content_width': 120},
    name='setup-analysis',
    help='Build a PyClone analysis configuration file.'
)
@click.option(
    '-i', '--in-files',
    multiple=True,
    required=True,
    type=click.Path(exists=True, resolve_path=True),
    help='''Path of TSV format file with copy number and allele count information. See build_mutations_file command
    for information. Maybe specified multiple times for multiple inputs.'''
)
@click.option(
    '-o', '--out-dir',
    required=True,
    type=click.Path(resolve_path=True),
    help='''Path of directory where config files will be placed.'''
)
@click.option(
    '-s', '--samples',
    multiple=True,
    help='''Sample name. Can be specified multiple times for multiple inputs. Should be in the same order as --in_files.
    If not set sample name will be inferred from file names and ordering in plots will be arbitrary.'''
)
@click.option(
    '-t', '--tumour-contents',
    multiple=True,
    type=float,
    help='''Tumour contents. Can be specified multiple times for multiple inputs. Should be in the same order as
    --in_files. If not set tumour content will be assumed to be 1.0.'''
)
@click.option(
    '-d', '--density',
    default='pyclone_beta_binomial',
    type=click.Choice(['pyclone_binomial', 'pyclone_beta_binomial']),
    help='''Tumour contents. Can be specified multiple times for multiple inputs. Should be in the same order as
    --in_files. If not set tumour content will be assumed to be 1.0.'''
)
@click.option(
    '--init-method',
    default='connected',
    type=click.Choice(['connected', 'disconnected']),
    help='''How to initialise the DP clustering algorithm. `connected` places all data points in one cluster, preferred
    for large datasets. `disconnected` places each data point in a separate cluster. Default `connected`.'''
)
@click.option(
    '--prior',
    default='major',
    type=click.Choice(['major', 'parental', 'total']),
    help='''Method used to set the possible genotypes. See online help for description. Default is major_copy_number.'''
)
@click.option(
    '--config-extras-file',
    default=None,
    type=click.Path(exists=True, resolve_path=True),
    help='''Path to YAML format file with additional or override config parameters. For advanced usage only.'''
)
def setup_analyis(**kwargs):
    run.setup_analysis(**kwargs)


#=======================================================================================================================
# Setup main interface
#=======================================================================================================================
@click.group()
def pyclone():
    pass


pyclone.add_command(build_table)
pyclone.add_command(plot_clusters)
pyclone.add_command(plot_loci)
pyclone.add_command(resume_analysis)
pyclone.add_command(run_analysis)
pyclone.add_command(setup_analyis)

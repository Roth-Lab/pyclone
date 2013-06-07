# License

PyClone is free for academic/non-profit use. For commercial use please contact sshah@bccrc.ca. Consult the LICENSE.txt 
file for more details. 

# Versions

## 0.12.1

* Fixed typos in some example files.

* Added command to plot parallel coordinates for multiple samples.

* Updated interface of plotting commands to take configuration files as arguments instead of traced directory.

## 0.12.0

* Changed input files to work from YAML config instead of command line arguments.

* Added ability to do multiple sample analysis.

* Added robust Beta-Binomial version of PyClone.

* Added genotype naive clustering DP methods with Gaussian, Binomial and Beta-Binomial densities.

* Updated and renamed the `build_inputs` -> `build_mutations_file` function for building YAML inputs from tsv file.

## 0.11.3

* Fixed overflow in mpear clustring.

## 0.11.1

* Small change to clustering to use mutation_id not mutation in output, to make consistent with simple input.

## 0.11

* Reverted to PyDP for implementing DP methods.

* Removed dependency on numpy in analysis code.

# Installation

To install PyClone make sure you have the necessary libraries (listed below) installed. After that PyClone installs like
any other Python package with `python setup.py install`.

If the installation worked correctly the `PyClone` command should now be available.

## Dependencies

### Required

The following packages are required to perform a basic analysis with PyClone.

* [PyDP >= 0.2.0](https://bitbucket.org/aroth85/pydp)

* [PyYAML >= 3.10](http://pyyaml.org)

### Optional

The following libraries are required to use the clustering and plotting capabilities of PyClone.

* [brewer2mpl >= 1.0] (https://github.com/jiffyclub/brewer2mpl) - Required for plotting.

* [eppl >= 0.1.0] (https://bitbucket.org/aroth85/eppl)

* [maplotlib >= 1.2.0](http://matplotlib.org) - Required for plotting.

* [numpy >= 1.6.2](http://www.numpy.org) - Required for plotting and clustering.

* [rpy2 >= 2.3.3](http://rpy.sourceforge.net/rpy2.html) - Only necessary to use the dynamic_tree_cut clustering method. The dynamicTreeCut tree cut package should also be installed in R.

* [scikits-learn >= 0.13](http://scikit-learn.org) - Only necessary to use affinity_propogation, dbscan, spectral_clustering clustering methods. 

* [scipy >= 0.11](http://www.scipy.org) - Required for plotting and clustering.

# Running PyClone

To run a PyClone analysis you need to perform several steps.

1. Prepare mutations input file(s).

2. Prepare a configuration file for the analysis.

3. Run the PyClone analysis using the `analyse` command.

4. (Optional) Plot results using the `plot_multi_sample`, `plot_cellular_frequencies` and `plot_similarity_matrix` commands.

5. (Optional) Cluster the PyClone output using the `cluster` command.

## Prepare An Input File

To run a PyClone analysis you need to prepare a properly formatted YAML file. An example file, "advance.yaml", ships with this package in the examples directory. To make it easier to produce such a file PyClone has a command `build_mutations_file` which will take a simpler tab separated file and produce the required YAML file.

### Simple Input

#### TSV Input File
The `build_mutations_file` takes a tab delimited file with a header as input and produces a YAML formatted file which can be used for running a PyClone analysis. An example input file 'simple.tsv' is included in the examples/ folder.

The required fields in this file are:

* mutation_id - A unique ID to identify the mutation. Good names are thing such a the genomic co-ordinates of the mutation i.e. chr22:12345. Gene names are not good IDs because one gene may have multiple mutations, in which case the ID is not unique and PyClone will fail to run.

* ref_counts - The number of reads covering the mutation which contain the reference (genome) allele.

* var_counts - The number of reads covering the mutation which contain the variant allele.

* normal_cn - The copy number of the cells in the normal population. For autosomal chromosomes this will be 2 and for sex chromosomes it could be either 1 or 2. For species besides human other values are possible.

* minor_cn - The minor copy number of the cancer cells. Usually this value will be predicted from WGSS or array data.

* major_cn - The major copy number of the cancer cells. Usually this value will be predicted from WGSS or array data.

NOTES: 

* If you do not major and minor copy number information you should set the minor copy number to 0, and the major copy number to the predicted total copy number. If you do this make sure to use the "total_copy_number" for the --var_prior flag of the `build_mutations_file` command. DO NOT use the "parental copy number" information method as it assumes you have knowledge of the minor and major copy number.

* Any additional columns in the tsv file will be ignored so feel free to add additional annotation fields.

#### From TSV -> YAML Mutations File

The information provided in this file is deliberately kept simple, but it is insufficient to run a PyClone analysis. In order to produce a file with enough information the `build_mutations_file` command attempts to guess some details about the possible states. All states guessed by this command will be weighted equally.

The key flags which provide some control over `build_mutations_file` are (typing PyClone build_mutations_file -h will also print out help):

* --ref_prior - Controls how the copy number of the cells in reference population is set. Setting the copy number of the reference population completely determines the genotype of the reference cells since all cells in the reference population have no variant alleles in their genotype by definition. The possibilities are to consider states where 1) "normal" - the reference population shares the same copy number as the normal population, 2) "variant" - there reference population shares the same copy number as the variant population, 3) "normal_variant" - states in which the reference population has the copy number of the normal or variant population are considered with equal probability. This only has an effect if the "no_zygosity" or "total_copy_number" methods are used with the --var_prior option.

* --var_prior - Method used to set the possible genotypes  of the variant population. "AB" assumes all mutations have the AB genotype. "BB" assumes all mutations have the BB genotype. "no_zygosity" assumes all mutation have the genotype with the predicted total copy number and one mutant allele i.e. AAB for a copy number 3 mutation. "parental_copy_number" sets considers all possible genotypes compatible with the predicted parental copy number. "total_copy_number" considers all possible genotypes compatible with the predicted total copy number. If reliable parental copy number is available the parental_copy_number method should be chosen. Default is total_copy_number.

### Advanced Input

For more control on specifying the states or prior weights a YAML file can be directly created. The file advanced.yaml in the examples/ directory shows the basic format. Rather than manually trying to format the file it is highly recommended that a YAML library such as PyYAML be used to help create the files. 

One top level nodes is required in the file.

* mutations - A list of mutations an there possible states. See below.

Under the mutations node we define each mutation as an item with the following nodes.

* id - The unique ID of the mutation.

* ref_counts - The number of reads covering the mutation which contain the reference (genome) allele.

* var_counts - The number of reads covering the mutation which contain the variant allele.

* states - A list of possible states for the sample at this mutation along with prior weights.

Under the states node we define the following elements.

* g_n - The genotype of the normal population in the state.

* g_r - The genotype of the reference population in the state.

* g_v - The genotype of the variant population in the state.

* prior_weight - The relative prior weight of the state. The values will be normalised across all states to create a valid probability.

## Building A PyClone Configuration File

Before you can run an PyClone analysis you need to prepare a YAML formatted configuration file. This file will specify the details of the MCMC run including the number of iterations, what model you want to use, and what samples need to be included.

The entries in the configuration file will vary depending on what value you set for the 'density' entry. The following entries will be in all configuration files. For example files see the examples/mixing directory.

### Common Entries

* working_dir - This is the working directory for the analysis. All paths in the configuration file will be relative to this. This is mainly use to reduce the verbosity of the configuration file.

* trace_dir - This is where the trace files for each parameter of the MCMC analysis will be written. This folder will be created if it does not exist.

* num_iters - The number of iterations of the MCMC chain.

* base_measure_params - This is a node which will have sub-entries depending on the density chosen. These sub-entries supply the parameters for the base measure of the Dirichlet process.

* concentration - This is a node which has several sub-entries related to the concentration parameter (alpha) for the Dirichlet Process.

	** value - Specifies the value for the concentration parameter. If this is to be inferred this value is simply the starting value.
	
	** priors - This node has two sub-entries. If this node is omitted the value specified for the concentration parameter will be fixed and not inferred.
	
		*** shape - The shape parameter in the Gamma prior on the concentration parameter.
		
		*** rate - The rate parameter in the Gamma prior on the concentration parameter.
		
* samples - This is a node which will contain one or more sub-entries which specify details about the samples used in the analysis. Each sub-entry is a node which contain a minimum of a unique sample ID and a path to TSV or YAML format input file. For the genotype naive densities, gaussian, binomial, and beta-binomial the input files follow the simple tsv format outline above. Copy number information will be ignored for these methods. For the pyclone_binomial or pyclone_beta_binomial methods you need to pass YAML formatted mutations files.

A sample input for genotype naive method is 

########################################################################################################################
samples:
  # Unique sample ID
  SRR385938:
    # Path where tsv formatted mutations file for the sample is placed.
    mutations_file: tsv/SRR385938.tsv
########################################################################################################################

and for a PyClone method is

########################################################################################################################
samples:
  # Unique sample ID
  SRR385938:
    # Path where YAML formatted mutations file for the sample is placed.
    mutations_file: yaml/parental_copy_number/SRR385938.yaml
    
    tumour_content:
      # The predicted tumour content for the sample. If you have no estimate set this to 1.0.
      value: 1.0
    
    # Expected sequencing error rate for sample
    error_rate: 0.001
########################################################################################################################      

## Run PyClone

Once the required YAML file has been created PyClone can be run using the `analyse` command. To see the possible flags for the command run `PyClone analyse -h`. The `analyse` command will run a full MCMC analysis of the data writing the results to several files in the specified output directory. All files are bz2 compressed tab separated files, so they can easily be read by external tools.

## Plot Results

After running the analyse you can use the `plot_multi_sample`, `plot_cellular_frequencies` and `plot_similarity_matrix` commands to visualise the results. These commands plot the prevalence (allelic or cellular) colour coded by cluster across samples, posterior density estimates of the cellular frequencies and a heatmap of the posterior similarity matrix.

## Cluster Results

PyClone provides some methods for producing flat clusterings of the posterior similarity matrix via the `cluster` command. The default method is "mpear" which is based on the clustering method defined in "Improved Criteria for Clustering Based on the Posterior Similarity Matrix" by Fritsch et al. Another method which may be of interest is the "dynamic_tree_cut" which was used in "The clonal and mutational evolution spectrum of primary triple-negative breast cancers" by Shah et al.
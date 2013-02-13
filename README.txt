# Installation

To install PyClone make sure you have the necessary libraries (listed below) installed. After that PyClone installs like
any other Python package with `python setup.py install`.

If the installation worked correctly the `PyClone` command should now be available.

## Dependencies

### Required

The following packages are required to perform a basic analysis with PyClone.

* [numpy 1.6.2](http://www.numpy.org)

* [PyYAML 3.10](http://pyyaml.org)

### Optional

The following libraries are required to use the clustering and plotting capabilities of PyClone.

* [brewer2mpl 1.0] (https://github.com/jiffyclub/brewer2mpl) - Required for plotting.

* [maplotlib 1.2.0](http://matplotlib.org) - Required for plotting.

* [rpy2 2.3.3](http://rpy.sourceforge.net/rpy2.html) - Only necessary to use the dynamic_tree_cut clustering method. The dynamicTreeCut tree cut package should also be installed in R.

* [scikits-learn 0.13](http://scikit-learn.org) - Only necessary to use affinity_propogation, dbscan, spectral_clustering clustering methods. 

* [scipy 0.11](http://www.scipy.org) - Required for plotting and clustering.

# Running PyClone

To run a PyClone analysis you need to perform several steps.

1. Prepare an input file.

2. Run the PyClone analysis using the `analyse` command.

3. (Optional) Plot results using the `plot_cellular_frequencies` and `plot_similarity_matrix` commands.

4. (Optional) Cluster the PyClone output using the `cluster` command.

## Prepare An Input File

To run a PyClone analysis you need to prepare a properly formatted YAML file. An example file, "advance.yaml", ships with this package in the examples directory. To make it easier to produce such a file PyClone has a command `build_input_file` which will take a simpler tab separated file and produce the required YAML file.

### Simple Input

The `build_input_file` takes a tab delimited file with a header as input and produces a YAML formatted file which can be used for running a PyClone analysis. An example input file 'simple.tsv' is included in the examples/ folder.

The required fields in this file are:

* mutation_id - A unique ID to identify the mutation. Good names are thing such a the genomic co-ordinates of the mutation i.e. chr22:12345. Gene names are not good IDs because one gene may have multiple mutations, in which case the ID is not unique.

* ref_counts - The number of reads covering the mutation which contain the reference (genome) allele.

* var_counts - The number of reads covering the mutation which contain the variant allele.

* cn_n - The copy number of the cells in the normal population. For autosomal chromosomes this will be 2 and for sex chromosomes it could be either 1 or 2. For species besides human other values are possible.

* cn_v - The copy number of the cells in the variant population. Usually this value will be predicted from WGSS or array data.

The information provided in this file is deliberately kept simple, but it is insufficient to run a PyClone analysis. In order to produce a file with enough information the `build_input_file` command attempts to guess some details about the possible states. All states guessed by this command will be weighted equally.

The key flags which provide some control over `build_input_file` are:

* --cn_r - Controls how the copy number of the cells in reference population is set. Setting the copy number of the reference population completely determines the genotype of the reference cells since all cells in the reference population have no variant alleles in their genotype by definition. The possibilities are to consider states where 1) "normal" - the reference population shares the same copy number as the normal population, 2) "variant" - there reference population shares the same copy number as the variant population, 3) "vague - states in which the reference population has the copy number of the normal or variant population are considered with equal probability.

* --g_v - Controls how to set the genotype for the variant population. Since we supply the copy number of the variant population in the input file, we only need to fill in information about how many variant alleles the genotype has. The possible choices are 1) "single" - assumes a single variant allele in the genotype i.e. AAAB, 2) "all" - assumes all alleles in the genotype are variant i.e. BBBB, 3) "vague" - considers all genotypes with a variant allele which are compatible with the predicted copy number assigning equal prior weight i.e. AAAB, AABB, ABBB, BBBB.  

### Advanced Input

For more control on specifying the states or prior weights a YAML file can be directly created. The file advanced.yaml in the examples/ directory shows the basic format. Rather than manually trying to format the file it is highly reccomemended that a YAML library such as PyYAML be used to help create the files. 

Two top level nodes are required in the file.

* error_rate - A scalar value indicating the sequencing error rate i.e. the probability of observing an A allele when it was really a B (assumed to be equal to the reverse).

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

## Run PyClone

Once the required YAML file has been created PyClone can be run using the `analyse` command. To see the possible flags for the command run `PyClone analyse -h`. The `analyse` command will run a full MCMC analysis of the data writing the results to several files in the specified output directory. All files are bz2 compressed tab separated files, so they can easily be read by external tools.

If the an estimate of tumour content is available be sure to set the --tumour_content flag.

## Plot Results

After running the analyse you can use the `plot_cellular_frequencies` and `plot_similarity_matrix` commands to visualise the results. These commands plot the posterior density estimates of the cellular frequencies and a heatmap of the posterior similarity matrix.

## Cluster Results

PyClone provides some methods for producing flat clusterings of the posterior similarity matrix via the `cluster` command. The default method is "mpear" which is based on the clustering method defined in "Improved Criteria for Clustering Based on the Posterior Similarity Matrix" by Fritsch et al. Another method which may be of interest is the "dynamic_tree_cut" which was used in "The clonal and mutational evolution spectrum of primary triple-negative breast cancers" by Shah et al.
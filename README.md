# License

PyClone is free for academic/non-profit use. For commercial use please contact sshah@bccrc.ca. Consult the LICENSE.txt 
file for more details. 

# Installation

## Using conda

You can install PyClone using bioconda.

`conda install pyclone -c bioconda -c conda-forge`

This will install PyClone into your current `conda` environment.
In some cases it may be better to create a separate `conda` environment for PyClone which be activated when needed.
This avoids issues due to conflicting libraries.
To create the environment execute the following command.

`conda create -n pyclone -c bioconda -c conda-forge pyclone`

Once the environment is created it can be activated using the following command.

`conda activate pyclone`

You can check that PyClone was installed correctly by running the following command which will show the help.

`PyClone --help`

## From source

PyClone is standard Python package.
You can find a list of dependencies in the `conda` recipe [here](https://github.com/Roth-Lab/pyclone/blob/master/conda/recipe/meta.yaml).
You will need to install PyDP from source as well, which can be found [here](https://github.com/Roth-Lab/pydp).

# Usage

## Input format

The majority of users will use PyClone by creating a set of tab delimited (tsv) input files, one file for each sample from the cancer.
The mandatory columns of this files are as follows.

- mutation_id: A unique identifier for the mutation. This should be the same across datasets.
- ref_counts: The number of reads overlapping the locus matching the reference allele.
- var_counts: The number of reads overlapping the locus matching the variant allele.
- normal_cn: The copy number of the locus in non-malignant cells. This should generally be 2 except for sex chromosomes in males.
- minor_cn: The copy number of the minor allele in the malignant cells. This must be less than equal the value in the major_cn column.
- major_cn: The copy number of the major allele in the malignant cells. This should be greater than equal to the value in the minor_cn column and greater than 0.

Any other columns will be ignored.
Example files are found [here](https://github.com/Roth-Lab/pyclone/tree/master/examples/mixing/tsv) from the mixing dataset used in the original PyClone paper. 

## Basic usage

The easiest way to run PyClone is using the `PyClone run_analysis_pipeline` pipeline command.
This will perform the steps to pre-process the input files, run the MCMC analysis and do the post-processing and plotting.
You will need to generate the input files as specified in the previous section.
You will need to pass two mandatory arguments.

- `--in_files`: A space delimited set of tsv files formatted as specified in the input format section.
- `--working_dir`: A directory where the pipeline will run and output results.

Two important optional flags are:

- `--tumour_contents`: A space delimited list of tumour content values between 0 and 1. The order of these should match the order samples were passed to `--in_files`. If this is not set tumour content is assumed to be 1. For most analysis it is important to set this or performance will suffer.
- `--samples`: A space delimited set of sample names with the order matching the order of `--in_files`. If this is not set the sample names will be inferred from the file names.

Additional arguments are available and can be listed using `PyClone run_analysis_pipeline --help`.

# Versions

## 0.13.1

- Added option to control how clusters are initialized with `--init_method` flag.
Using `--init_method connected` will place all data points in the same cluster leading to faster sampling at the risk of getting stuck in local modes.

- Added option to control the maximum number of clusters considered during post-processing using `--max_clusters` flag.
For large datasets this should be set to values less than 100 to speed things up.

- Using [numba](https://numba.pydata.org/) to speed up code.

## 0.13.0

Most changes in this release are internal refactoring of the code and should be invisble to the user.

- Removed IGMM, IBMM, and IBBMM methods

- Changed plotting code to use seaborn
	- Removed eppl dependency
	- Removed brewer2mpl dependency

- Changed output of build_table to a tidy data frame
	- Old style table available with the --old_style flag
	
- Re-wrote the multi-sample plotting code
	- Should fix issues with too many clusters
	- Slightly nicer looking plots

- Switched to using pandas for data wrangling

- Changed nomenclature.
	- Variant allele frequency (VAF): proportion of reads with variant
	- Cellular prevalence: proportion of cancer cells harbouring a mutation
	- The name of the trace files have been altered
	- The name of some function calls have been altered
	- Plot labels have been altered

## 0.12.9

* Fixed a bug in the calculation of the mpear code.

## 0.12.8

* Fixed a bug causing PyClone BetaBinomial to not support fixed precision parameter.

## 0.12.7

* Fixed bug in multi-sample plotting of allelic prevalences. Note this requires an upgrade to eppl 0.2.3 or greater.

## 0.12.6

* Modified plot_cellular_frequencies to accept the `--file_format` flag which sets the output format for plot files.

## 0.12.5

* Changed multi sample plot to show error bars

* The `build_table` no outputs the std error of the cellular prevalence trace for each mutation in each sample 

## 0.12.4

* Fixed bug in the multi-sample plotting/table code which caused failures for non PyClone densities.

## 0.12.3

* Added ability to output table used for generating multi-sample plots.

## 0.12.2

* Fixed a bug which would cause `build_mutations_file` to fail if the output was given as a relative path.

* Updated the configuration files in the examples/ folder to have saner values.

* Changed README to redirect to website.

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

* Fixed overflow in mpear clustering.

## 0.11.1

* Small change to clustering to use mutation_id not mutation in output, to make consistent with simple input.

## 0.11

* Reverted to PyDP for implementing DP methods.

* Removed dependency on numpy in analysis code.

## Older

* Unfortunately I did not keep a complete list of changes as the software evolved.

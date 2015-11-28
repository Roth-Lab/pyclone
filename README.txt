# License

PyClone is free for academic/non-profit use. For commercial use please contact sshah@bccrc.ca. Consult the LICENSE.txt 
file for more details. 

Please visit https://bitbucket.org/aroth85/pyclone for installation and usage help.

# Versions

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

## 0.12.10

* Internal fix to use mpear implementation in PyDP. Requires PyDP >= 0.2.3

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
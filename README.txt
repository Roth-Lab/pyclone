# License

PyClone is free for academic/non-profit use. For commercial use please contact sshah@bccrc.ca. Consult the LICENSE.txt 
file for more details. 

Please visit https://bitbucket.org/aroth85/pyclone for installation and usage help.

# Versions

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
# Installation

To install PyClone make sure you have the necessary libraries (listed below) installed. After that PyClone installs like
any other Python package with `python setup.py install`.

If the installation worked correctly the `PyClone` command should now be available.

## Dependencies

### Required

The following packages are required to perform a basic analysis with PyClone.

* [numpy 1.6.2](http://www.numpy.org)


### Optional

The following libraries are required to use the plotting capabilities of PyClone.

* [brewer2mpl 1.0] (https://github.com/jiffyclub/brewer2mpl)

* [maplotlib 1.2.0](http://matplotlib.org)

* [rpy2 2.3.3](http://rpy.sourceforge.net/rpy2.html)

* [scikits-learn 0.13](http://scikit-learn.org)

* [scipy 0.11](http://www.scipy.org)

# Running PyClone

To run a PyClone analysis you need to perform several steps.

1. Prepare an input file.

2. Run the PyClone analysis using the `analyse` command.

3. (Optional) Plot results using the `plot_cellular_frequencies` and `plot_similarity_matrix` commands.

4. (Optional) Cluster the PyClone output using the `cluster` command.

## Prepare an input file

The input file for PyClone is a tab separated file with a header. The header should contain the following columns in any
order. Note any additional columns will be *ignored* by PyClone so its save to include them to help document the file.

* mutation - This column has the unique ID for the mutation. Using gene names is a bad idea since a single gene may have
multiple mutations. If the same ID is used twice PyClone will raise an error.

* b - The number of reads covering the mutant loci which contain the variant allele.

* d - The total number of reads covering the mutant loci.

* cn_r - A comma separated list of copy numbers for the reference population. This should be the same length as the list
for cn_v, mu_v and prior_weight.

* cn_v -  A comma separated list of copy numbers for the variant population. This should be the same length as the list
for cn_r, mu_v and prior_weight.

* mu_v -  A comma separated list of probabilities from sampling a variant allele from a cell with a given genotype. For 
example if the possible genotypes are AB and BB the entry would be 0.5, 0.999 (assumin an error rate of 0.001. This 
should be the same length as the list for cn_r, cn_v and prior_weight.

* prior_weight - A comma separated list of relative prior weights assigned to a given state of the sample. The entries 
in this last match with the entries in cn_r, cn_v, mu_v so it must be the same length.

An example file `pyclone.example.tsv` is included under the examples directory.

To understand the input consider the following lines

mutation	b	d	cn_r	cn_v	mu_v	prior_weight
mutation_2	500	1000	2, 2	2,3	0.999,0.666	1,2

The first line is the header of the file. The second line says we have a mutation called `mutation_2`. There are b=500
reads containing the variant allele covering the loci. There are d=1000 total reads covering the loci.

Now notice that cn_r, cn_v, mu_v and prior_weight are all lists of length 2. The first entry of each list go together to 
form a state and its prior weight. Similarly the second entries go together to form a state and ist prior weight.

The first state assumes: 

* The reference population has copy number 2. 

* The variant population has copy number 2.

* The probability of sampling a variant read from the cells in the variant population is 0.999. This would be compatible
with a genotype of BB assuming there was an error rate of 0.001 for sequencing.

* The prior weight for this state is 1.

The second state assumes:

* The reference population has copy number 2.

* The variant population has copy number 3.

* The probability of sampling a variant read from the cells in the variant population is 0.666~2/3. This would be 
compatible with a genotype of ABB. Note there are 2 B alleles out of 3 total alleles in the genotype.

* The prior weight for this state is 2.

## Prior Weights

The prior weights represent our belief ahead of time about the state of the sample. These weights are relative and will
be normalised by PyClone to sum to 1 to form valid prior probabilities.


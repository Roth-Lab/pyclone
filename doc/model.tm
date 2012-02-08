<TeXmacs|1.0.7.14>

<style|generic>

<\body>
  <doc-data|<doc-title|PyClone : Software For Inferring Cellular Frequencies
  From Allellic Count Data>|<doc-author-data|<author-name|Andrew Roth>>>

  <section|Introduction>

  <subsection|Notation>

  In what follows we assume that we have aligned sequence data. For each
  position of interest this data is summarised by counting how many reads
  contain a given nucleotide from the set
  <math|\<Sigma\>=<around*|{|A,C,T,G|}>>. For simplicity we further assume
  that there is a reference genome used and that there are only two alleles
  at a given position, the reference allele
  <math|<with|mode|text|A>\<in\>\<Sigma\>> and the variant allele
  <math|<with|mode|text|B>\<in\>\<Sigma\>>. For each position, <math|i>, we
  reduce read data to count data, <math|a<rsup|i>>, the number of reads with
  nucleotides matching this reference allele and <math|b<rsup|i>> the number
  of reads with nucleotides matching the variant allele.

  With this formalism a diploid genome can have one of three possible
  genotypes, <math|\<nocomma\>\<nocomma\><with|mode|text|AA>,<with|mode|text|AB>\<nocomma\>,<with|mode|text|BB>>,
  at each position. Cancer genomes are not strictly diploid, so we need to
  consider an extended set of possible genotypes. For each position the
  genotype can be fully specified by two value. The number of copies of that
  position, <math|c<rsup|i>>, and the number reference allelels in the
  genotype <math|r<rsup|i>>. For example the genotype
  <math|<with|mode|text|AAAB>> would have <math|c<rsup|i>=4> and
  <math|r<rsup|i>=3>.

  In the absence of any other confounding factors, the frequency of reference
  alleles we observe at a given position should be
  <math|f<rsup|i>=<frac|r<rsup|i>|c<rsup|i>>>. Due to sampling and sequencing
  error the observed reference alleleic frequency <math|<wide|f|^><rsup|i>>
  will generally not be exactly equal to <math|f<rsup|i>>. However, as depth
  of coverage increases <math|<wide|f|^><rsup|i>> should converge to
  <math|f<rsup|i>>.

  <subsection|Sequencing Error>

  Sequencing error means that in general <math|<wide|f|^><rsup|i>> will
  differ slightly from <math|f<rsup|i>>. If the errors are random, then for
  genotypes where <math|r<rsup|i>\<neq\>0> or
  <math|r<rsup|i>\<neq\>c<rsup|i>> the errors should cancel out and this
  effect should not be noticeable. However, in the two exceptional cases
  where <math|r<rsup|i>=0> or <math|r<rsup|i>=c<rsup|i>> we assume there will
  be a small random deviation, <math|\<varepsilon\>\<ll\>1>, in observed
  frequency assosciated with the error rate of the sequencing technology.
  That is when <math|r<rsup|i>=0> we should expect
  <math|f<rsup|i>=\<varepsilon\>> and when <math|r<rsup|i>=c<rsup|i>> then
  <math|f<rsup|i>=1-\<varepsilon\>>.

  <subsection|Heterogeneity>

  We define heterogeneity to mean there are two or more populations of cells
  with different genotypes at position <math|i>. Heterogeneity at a locus
  means that <math|f<rsup|i>> is no longer well defined, the
  <math|<wide|f|^><rsup|i>> will in general not directly inform us about the
  underlying genotypes.\ 

  For example we can imagine a simple example where <math|50%> of the cells
  have the genotype <math|<with|mode|text|AA>> and <math|50%> have the
  genotype <math|<with|mode|text|BB>>. Then at high depth, we expect
  <math|<wide|f|^><rsup|i>> to converge to
  <math|0.5\<cdot\>f<rsub|<with|mode|text|AA>>+0.5\<cdot\>f<rsub|<with|mode|text|AB>>>,
  where with some abuse of notation we use <math|f<rsub|g>> to be the allelic
  frequency associated with genotype <math|g>.

  If we assume there are <math|J> sub-populations at site <math|i> we can
  label the genotypes as using <math|c<rsup|i><rsub|j>,r<rsub|j><rsup|i>,f<rsub|j><rsup|i>>
  to indicate the copy number, number of reference alleles and associated
  reference allele frequency in the <math|j<rsup|th>> sub-population. Note
  the if we have <math|r<rsup|i><rsub|j>=0> or
  <math|r<rsub|j><rsup|i>=c<rsub|j><rsup|i>>, for all <math|j>, the
  heterogeneity becomes unobservable since <math|f<rsub|j><rsup|i>> will be
  <math|\<varepsilon\>> or <math|1-\<varepsilon\>> for all populations.

  Again an example may clarify this. Consider a position with two
  sub-populations each at <math|50%> with genotypes
  <math|<with|mode|text|AA>> and <math|<with|mode|text|AAA>>. Then
  <math|f<rsub|<with|mode|text|AA>>=f<rsub|<with|mode|text|AAA>>=1-\<varepsilon\>>,
  so that <math|<wide|f|^><rsup|i>> will converge to
  <math|0.5\<cdot\>f<rsub|<with|mode|text|AA>>+0.5\<cdot\>f<rsub|<with|mode|text|AAA>>=0.5\<cdot\>
  <around*|(|1-\<varepsilon\>|)>+0.5\<cdot\><around*|(|1-\<varepsilon\>|)>=1-\<varepsilon\>>.

  This observation is important because it means that allelic frequencies are
  not useful for distinguishing copy number variants.

  <section|Model>

  <subsection|Assumptions>

  To simplify the subsequent analysis we make a simplifying assumption that
  there are two populations of cells at a given position <math|i>. One
  population consits of cells for which <math|r<rsup|i><rsub|j>=c<rsub|j><rsup|i>>,
  which we refer to as the <em|reference> population. This includes both
  normal cells, and tumour cells with copy number variations. The common
  factor in this group is that the genotypes contain no variant alleles. As
  discussed above, using allele frequencies we cannot deconvolute the
  fraction of normal cells, and tumour cells with no mutations. We assign the
  variable <math|f<rsub|r><rsup|i>=1-\<varepsilon\>> to the reference allele
  frequency for this population.

  The second population which we refer to as the <em|variant> population will
  have <math|r<rsub|j><rsup|i>\<less\>c<rsub|j><rsup|i>>, that is at least
  one variant allele is present in the genotype. We will further assume that
  only one such population exists, so that we have the parameters
  <math|c<rsub|v><rsup|i>,r<rsub|v><rsup|i>,f<rsub|v><rsup|i>> for this
  population.

  We let <math|\<phi\><rsup|i>> be the \ fraction of cells from the variant
  population and <math|1-\<phi\><rsup|i>> the fraction of cells from the
  reference population. It then follows that <math|<wide|f|^><rsup|i>>
  converges to <math|<around*|(|1-\<phi\><rsup|i>|)>f<rsub|r><rsup|i>+\<phi\><rsup|i>
  f<rsub|v><rsup|i>>. This implies that if we knew <math|f<rsub|v><rsup|i>>,
  we could estimate <math|\<phi\><rsup|i>> via

  <\eqnarray>
    <tformat|<table|<row|<cell|\<phi\><rsup|i>>|<cell|=>|<cell|<frac|<wide|f|^><rsup|i>-f<rsub|r><rsup|i>|f<rsub|v><rsup|i>-f<rsub|r><rsup|i>>>>>>
  </eqnarray>

  The key difficulty is that knowledge of <math|f<rsub|v><rsup|i>>, would
  imply knowledge of the variant populations genotype. This information is
  generally not available, however it is possible that some prior beliefs
  over possible genotypes exists.

  <subsection|Accomodating Genotype Uncertainty>

  The discussion above implies that we can model the number of reference
  allele counts observed at position <math|i>, as a binomial distribution
  with paramters <math|n=d<rsup|i>=a<rsup|i>+b<rsup|i>>, and
  <math|p=<around*|(|1-\<phi\><rsup|i>|)>f<rsub|r><rsup|i>+\<phi\><rsup|i>
  f<rsub|v><rsup|i>>. However, the unceratainty in genotype means in
  <math|f<rsub|v><rsup|i>> is unknown.\ 

  At this point it is useful to switch notation slightly. We will let
  <math|\<mu\><rsub|r>=1-\<varepsilon\>> be the probability of sampling a
  reference allele from the reference population. To reiteratre,
  <math|\<varepsilon\>\<ll\>1> is a value assosciated with the error rate of
  the sequencing technology.

  We will also introduce an infinite vector
  <math|\<b-mu\><rsub|v>=<around*|(|\<mu\><rsub|<with|mode|text|A>>,\<mu\><rsub|<with|mode|text|B>>,\<mu\><rsub|<with|mode|text|AA>>,\<ldots\>|)>>.
  The entries, <math|\<mu\><rsub|v:g>>, in this vector are the probability of
  sampling a reference allele from the variant population with genotype
  <math|g>.\ 

  In addition for each position <math|i> we will introduce a vectors
  <math|\<b-pi\><rsup|i>=<around*|(|\<pi\><rsup|i><rsub|<with|mode|text|A>>,\<pi\><rsup|i><rsub|<with|mode|text|B>>,\<pi\><rsup|i><rsub|<with|mode|text|AA>>,\<ldots\>|)>>,
  in which the entries, <math|\<pi\><rsup|i><rsub|g>>, are the probabilities
  the variant population at site <math|i> has genotype <math|g>. Only a
  finite number of entries in <math|\<b-pi\><rsup|i>> will be non-zero.

  Finally we introduce a variable <math|G<rsup|i>>, which indicates which
  genotype the variant population has.

  We can the create a hierachical bayesian model of the data as follows.

  <\eqnarray>
    <tformat|<table|<row|<cell|\<phi\><rsup|i>>|<cell|\<sim\>>|<cell|<with|mode|text|Uniform><around*|(|0,1|)>>>|<row|<cell|G<rsup|i>>|<cell|\<sim\>>|<cell|<with|mode|text|Discrete><around*|(|\<b-pi\><rsup|i>|)>>>|<row|<cell|a<rsup|i>\|d<rsup|i>,G<rsup|i>=g,\<phi\><rsup|i>,\<mu\><rsub|r>,\<b-mu\><rsub|v>>|<cell|\<sim\>>|<cell|<with|mode|text|Binomial><around*|(|d<rsup|i>,<around*|(|1-\<phi\><rsup|i>|)>\<mu\><rsub|r>+\<phi\><rsup|i>
    \<mu\><rsub|v:g>|)>>>>>
  </eqnarray>

  We can then compute a posterior distribution on <math|\<phi\><rsup|i>>

  <\eqnarray>
    <tformat|<table|<row|<cell|\<bbb-P\><around*|(|\<phi\><rsup|i>,G<rsup|i>=g\|a<rsup|i>,d<rsup|i>,\<phi\><rsup|i>,\<mu\><rsub|r>,\<b-mu\><rsub|v>,\<b-pi\><rsup|i>|)>>|<cell|\<propto\>>|<cell|\<bbb-P\><around*|(|a<rsup|i>\|d<rsup|i>,G<rsup|i>=g,\<phi\><rsup|i>,\<mu\><rsub|r>,\<b-mu\><rsub|v>|)>\<bbb-P\><around*|(|G<rsup|i>\|\<b-pi\><rsup|i>|)>\<bbb-P\><around*|(|\<phi\><rsup|i>|)>>>>>
  </eqnarray>

  We can then marginalise out <math|G<rsup|i>> to obtain a posterior
  distribution over <math|\<phi\><rsup|i>>.

  <subsection|Hierachical Modelling Of Genotype Uncertainty>

  Eliciting prior information about <math|\<b-pi\><rsup|i>> will be
  difficult. It is beneficial to add another level to the hierachy of the
  model by placing a Dirichle prior over <math|\<b-pi\><rsup|i>>. This frees
  us from having to directly input prior probabilities about genotypes, and
  instead allows us to work with the simpler pseudo-count parameters of the
  Dirichlet distribution. In addition, depeening the hierachy will provide
  some measure of protection against inaccurate prior specification for
  <math|\<b-pi\>>.

  Formally we have

  <\eqnarray>
    <tformat|<table|<row|<cell|\<b-pi\><rsup|i>\|\<b-delta\><rsup|i>>|<cell|\<sim\>>|<cell|<with|mode|text|Dirichlet><around*|(|\<b-delta\><rsup|i>|)>>>>>
  </eqnarray>

  Calculation of the posterior then becomes

  <\eqnarray>
    <tformat|<table|<row|<cell|\<bbb-P\><around*|(|\<phi\><rsup|i>,G<rsup|i>=g,\<b-pi\><rsup|i>\|a<rsup|i>,d<rsup|i>,\<phi\><rsup|i>,\<mu\><rsub|r>,\<b-mu\><rsub|v>,\<b-delta\><rsup|i>|)>>|<cell|\<propto\>>|<cell|\<bbb-P\><around*|(|a<rsup|i>\|d<rsup|i>,G<rsup|i>=g,\<phi\><rsup|i>,\<mu\><rsub|r>,\<b-mu\><rsub|v>|)>\<bbb-P\><around*|(|G<rsup|i>\|\<b-pi\><rsup|i>|)>\<bbb-P\><around*|(|\<phi\><rsup|i>|)>\<bbb-P\><around*|(|\<b-pi\><rsup|i>\|\<b-delta\><rsup|i>|)>>>>>
  </eqnarray>

  Again we can margianlise the nuisance parameter <math|G<rsup|i>>, and now
  <math|\<b-pi\><rsup|i>>. The final form of the posterior after
  mariganalisation is

  <\eqnarray>
    <tformat|<table|<row|<cell|\<bbb-P\><around*|(|\<phi\><rsup|i>\|a<rsup|i>,d<rsup|i>,\<phi\><rsup|i>,\<mu\><rsub|r>,\<b-mu\><rsub|v>,\<b-delta\><rsup|i>|)>>|<cell|\<propto\>>|<cell|<big|sum><rsub|g><around*|[|<frac|<big|prod><rsub|g<rprime|'>\<neq\>g>
    \<Gamma\><around*|(|\<delta\><rsub|g<rprime|'>>|)>\<times\>\<Gamma\><around*|(|\<delta\><rsub|g>+1|)>|\<Gamma\><around*|(|<big|sum><rsub|g<rprime|'>>\<delta\><rsub|g<rprime|'>>+1|)>>|]><with|mode|text|Binomial><around*|(|d<rsup|i>,<around*|(|1-\<phi\><rsup|i>|)>\<mu\><rsub|r>+\<phi\><rsup|i>
    \<mu\><rsub|v:g>|)>>>>>
  </eqnarray>

  Since <math|\<phi\><rsup|i>> is a one dimensional parameter with support in
  the interval <math|<around*|[|0,1|]>> it is trivial to compute the
  normalisation constant using numerical integration techniques.

  <subsection|Sharing Statistical Strength Across Samples>

  In the above formulation piror uncertainty about variant population
  genotypes, will be translated into posterior uncertainty about
  <math|\<phi\><rsup|i>>. In particular, if we believe several genotypes are
  equally probable, then there will be equal number of modes of the same
  heigh in the posterior for <math|\<phi\><rsup|i>> corresponding to each
  genotype. Without further assumptions we have no reason to believe any of
  these modes are more probable than the others.

  If we consider the set of all mutations in a sample together, this no
  longer need be true. By utilising a shared prior for the cellular
  frequencies <math|\<phi\><rsup|i>>, across all positions in the sample we
  can impose extra constraints to resolve ambiguities in genotypes. The
  choice of such a prior distribution will have a dramatic effect on
  inference, and should be considered carefully. In particular, the
  biological implications of the prior must be fully understood.

  <subsubsection|Dirichlet Process Prior>

  If we are willing to assume that modes the posterior of
  <math|\<phi\><rsup|i>> which are common between positions are more likely,
  we can resolve the multimodal posteriors. Biologically these means we are
  assuming there are groups of mutations, which tend to occur at the same
  frequency.\ 

  We can explictily incorporate such an assumption into the model, if we
  change the prior distribution of <math|\<phi\><rsup|i>> to be shared across
  all positions in a sample. In particular we can use a distribution over a
  discrete set of atom for all <math|\<phi\><rsup|i>> in the sample. Using
  simple discrete distribution, and adjusting the location of the atoms would
  entail selecting the number of groups of mutation in advance. Again, this
  information is not generally available. We can avoid this difficulty by
  using a semi-parameteric Dirichlet process prior (DPP) over
  <math|\<phi\><rsup|i>>.

  \;

  \;

  \;
</body>

<\initial>
  <\collection>
    <associate|page-type|letter>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-10|<tuple|2.4.1|?>>
    <associate|auto-2|<tuple|1.1|?>>
    <associate|auto-3|<tuple|1.2|?>>
    <associate|auto-4|<tuple|1.3|?>>
    <associate|auto-5|<tuple|2|?>>
    <associate|auto-6|<tuple|2.1|?>>
    <associate|auto-7|<tuple|2.2|?>>
    <associate|auto-8|<tuple|2.3|?>>
    <associate|auto-9|<tuple|2.4|?>>
    <associate|footnote-1|<tuple|1|?>>
    <associate|footnr-1|<tuple|1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Introduction>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1.5fn>|Notation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|3fn>|Sequencing Error
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|3fn>|Heterogeneity
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|3fn>|Assumptions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|3fn>|Accomodating Genotype Uncertainty
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>
    </associate>
  </collection>
</auxiliary>
# generateMhlReport

This function computes *Linearised* Methylated Haplotype Load (\\lMHL\\)
per genomic position.

## Usage

``` r
generateMhlReport(
  bam,
  report.file = NULL,
  haplotype.context = c("CG", "CHG", "CHH", "CxG", "CX"),
  max.haplotype.window = 0,
  min.haplotype.length = 0,
  max.outofcontext.beta = 0.1,
  ...,
  gzip = FALSE,
  verbose = TRUE
)
```

## Arguments

- bam:

  BAM file location string OR preprocessed output of
  [`preprocessBam`](preprocessBam.md) function. Read more about BAM file
  requirements and BAM preprocessing at
  [`preprocessBam`](preprocessBam.md).

- report.file:

  file location string to write the \\lMHL\\ report. If NULL (the
  default) then report is returned as a
  [`data.table`](https://rdatatable.gitlab.io/data.table/reference/data.table.html)
  object.

- haplotype.context:

  string for a cytosine context that defines a haplotype:

  - "CG" (the default) – CpG cytosines only (called as zZ)

  - "CHG" – CHG cytosines only (xX)

  - "CHH" – CHH cytosines only (hH)

  - "CxG" – CG and CHG cytosines (zZxX)

  - "CX" – all cytosines; this, as well as the other non-CG contexts,
    may have little sense but still included for consistency

  If \\lMHL\\ calculations are needed for all three possible cytosine
  contexts *independently*, one has to run this function for each
  required \`haplotype.context\` separately, because
  \`haplotype.context\`=="CX" assumes that *any* cytosine context is
  allowed within the same haplotype. This behaviour may change in the
  future.

- max.haplotype.window:

  non-negative integer for maximum value of \\L'\\ in \\lMHL\\ formula.
  When 0 (the default), calculations are performed for the full
  haplotype length (\\L'=L\\, although the maximum value is currently
  limited to 65535). Having no length restrictions make sense for
  short-read sequencing when the length of the read is comparable to the
  length of a typical methylated block, the depth of coverage is high,
  and the lengths of all reads are roughly equal. However, calculations
  using non-restricted haplotype length are meaningless for long-read
  sequencing — when the same read may cover a number of regions with
  very different methylation properties, and reads themselves can be of
  a very different length. In the latter case it is advised to limit the
  \`max.haplotype.window\` to a number of cytosines in a typical
  hypermethylated region. For thorough explanation and more examples,
  see Details section and vignette.

- min.haplotype.length:

  non-negative integer for minimum length of a haplotype (default: 0
  will include haplotypes of any length). When
  \`min.haplotype.length\`\>0, reads (read pairs) with fewer than
  \`min.haplotype.length\` cytosines within the \`haplotype.context\`
  are skipped.

- max.outofcontext.beta:

  real number in the range \[0;1\] (default: 0.1). Reads (read pairs)
  with average beta value for out-of-context cytosines **above** this
  threshold are skipped. Set to 1 to disable filtering.

- ...:

  other parameters to pass to the [`preprocessBam`](preprocessBam.md)
  function. Options have no effect if preprocessed BAM data was supplied
  as an input.

- gzip:

  boolean to compress the report (default: FALSE).

- verbose:

  boolean to report progress and timings (default: TRUE).

## Value

[`data.table`](https://rdatatable.gitlab.io/data.table/reference/data.table.html)
object containing \\lMHL\\ report or NULL if report.file was specified.
The report columns are:

- rname – reference sequence name (as in BAM)

- strand – strand

- pos – cytosine position

- context – methylation context

- coverage – number of reads (read pairs) that include this position

- length – average length of a haplotype, i.e., average number of
  cytosines within \`haplotype.context\` for reads (read pairs) that
  include this position

- lmhl – \\lMHL\\ value

## Details

The function reports *Linearised* Methylated Haplotype Load (\\lMHL\\)
at the level of individual cytosines using BAM file location or
preprocessed data as an input. Function uses the following formula:

\$\$lMHL=\frac{\sum\_{i=1}^{L'} w\_{i} \times MH\_{i}}{\sum\_{i=1}^{L'}
w\_{i} \times H\_{i}}\$\$

where \\L'\\ is the length of a calculation window (e.g., number of
CpGs; \\L' \le L\\, where \\L\\ is the length of a haplotype covering
current genomic position), \\MH\_{i}\\ is a number of fully successive
methylated stretches with \\i\\ loci within a methylated stretch that
overlaps current genomic position, \\H\_{i}\\ is a number of fully
successive stretches with \\i\\ loci, \\w\_{i}\\ is a weight for
\\i\\-locus haplotype (\\w\_{i}=i\\).

This formula is a modification of the original Methylated Haplotype Load
(MHL) formula that was first described by Guo et al., 2017 (doi:
[10.1038/ng.3805](https://doi.org/10.1038/ng.3805)):

\$\$MHL=\frac{\sum\_{i=1}^{L} w\_{i} \times
\frac{MH\_{i}}{H\_{i}}}{\sum\_{i=1}^{L} w\_{i}}\$\$

where \\L\\ is the length of a longest haplotype covering current
genomic position, \\\frac{MH\_{i}}{H\_{i}}=P(MH\_{i})\\ is the fraction
of fully successive methylated stretches with \\i\\ loci, \\w\_{i}\\ is
a weight for \\i\\-locus haplotype (\\w\_{i}=i\\).

The modifications to original formula are made in order to:

- **provide granularity of values** — the original MHL formula gives the
  same MHL value for every cytosine of a partially methylated haplotype
  (e.g., MHL=0.358 for each cytosine within a read with methylation call
  string "zZZZ"). In contrast, \\lMHL\\==0 for the non-methylated
  cytosines (e.g., \\lMHL\\==c(0, 0.5, 0.5, 0.5) for cytosines within a
  read with methylation call string "zZZZ").

- **enable calculations for long-read sequencing alignments** — \\lMHL\\
  calculation window can be limited to a particular number of cytosines.
  This allows to use the formula for very long haplotypes as well as to
  compare values for sequencing data of varying read length.

- **reduce the complexity of MHL calculation** for data of high breadth
  and depth — \\lMHL\\ values for all genomic positions can be
  calculated using a single pass (cycling through reads just once) as
  the linearised calculations of numerator and denominator for \\lMHL\\
  do not require prior knowledge on how many reads cover a particular
  position. This is achieved by moving \\H\_{i}\\ multiplier to the
  denominator of the \\lMHL\\ formula.

These modifications make \\lMHL\\ calculation similar though
*non-equivalent* to the original MHL. However, the most important
property of MHL — emphasis on hypermethylated blocks — is retained. And
in return, \\lMHL\\ gets better applicability for analysis of sequencing
data of varying depth and read length.

Other notes on function's behaviour:

Methylation string bases in unknown context ("uU") are simply ignored,
which, to the best of our knowledge, is consistent with the behaviour of
other tools.

Cytosine context present in more than 50% of the reads is assumed to be
correct, while all bases at the same position but having other
methylation context are simply ignored. This allows reports to be
prepared without using the reference genome sequence.

## See also

\`values\` vignette for a comparison and visualisation of epialleleR
output values for various input files. \`epialleleR\` vignette for the
description of usage and sample data.

[`preprocessBam`](preprocessBam.md) for preloading BAM data,
[`generateCytosineReport`](generateCytosineReport.md) for other
methylation statistics at the level of individual cytosines,
[`generateBedReport`](generateBedReport.md) for genomic region-based
statistics, [`generateVcfReport`](generateVcfReport.md) for evaluating
epiallele-SNV associations, [`extractPatterns`](extractPatterns.md) for
exploring methylation patterns and [`plotPatterns`](plotPatterns.md) for
pretty plotting of its output, [`generateBedEcdf`](generateBedEcdf.md)
for analysing the distribution of per-read beta values.

## Examples

``` r
  capture.bam <- system.file("extdata", "capture.bam", package="epialleleR")
  
  # lMHL report
  mhl.report <- generateMhlReport(capture.bam)
#> Checking BAM file: 
#> short-read, paired-end, name-sorted alignment detected
#> Reading paired-end BAM file 
#> [0.012s]
#> Preparing lMHL report 
#> [0.021s]
  
  # lMHL report with a `max.haplotype.window` of 1 is identical to a
  # conventional cytosine report (or nearly identical when sequencing errors
  # are present)
  mhl.report <- generateMhlReport(capture.bam, max.haplotype.window=1)
#> Checking BAM file: 
#> short-read, paired-end, name-sorted alignment detected
#> Reading paired-end BAM file 
#> [0.012s]
#> Preparing lMHL report 
#> [0.020s]
  cg.report  <- generateCytosineReport(capture.bam, threshold.reads=FALSE)
#> Checking BAM file: 
#> short-read, paired-end, name-sorted alignment detected
#> Reading paired-end BAM file 
#> [0.012s]
#> Preparing cytosine report 
#> [0.012s]
  identical(
    mhl.report[, .(rname, strand, pos, context, value=lmhl)],
    cg.report[ , .(rname, strand, pos, context, value=meth/(meth+unmeth))]
  )
#> [1] TRUE
```

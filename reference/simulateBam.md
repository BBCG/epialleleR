# simulateBam

This function creates sample BAM files given mandatory and optional BAM
fields.

## Usage

``` r
simulateBam(
  output.bam.file = NULL,
  qname = NULL,
  flag = NULL,
  rname = NULL,
  pos = NULL,
  mapq = NULL,
  cigar = NULL,
  rnext = NULL,
  pnext = NULL,
  tlen = NULL,
  seq = NULL,
  qual = NULL,
  ...,
  verbose = TRUE
)
```

## Arguments

- output.bam.file:

  output BAM file location string. If NULL (default), records are not
  written to BAM but returned as a
  [`data.table`](https://rdatatable.gitlab.io/data.table/reference/data.table.html)
  object for review.

- qname:

  character vector of query names. When default (NULL), names like
  "q0001".."qNNNN" will be assigned.

- flag:

  integer vector of bitwise flags (a combination of the BAM_F\*
  constants). When default (NULL), zero (i.e., unique, valid,
  single-end, aligned read) is assigned for every record.

- rname:

  character vector of chromosome (reference) names. When default (NULL),
  "chrS" is assigned for every record.

- pos:

  integer vector of 1-based leftmost coordinates of the queries. When
  default (NULL), 1 is assigned for every record.

- mapq:

  integer vector of mapping qualities. When default (NULL), 60 is
  assigned for every record.

- cigar:

  character vector of CIGAR strings. When default (NULL), "lM" is
  assigned for every record, where \`l\` is the length of the query
  (\`seq\`).

- rnext:

  character vector of chromosome (reference) names for next read in
  template. When default (NULL), "chrS" is assigned for every record.

- pnext:

  integer vector of 1-based leftmost coordinates of next read in
  template. When default (NULL), 1 is assigned for every record.

- tlen:

  integer vector of observed template lengths. When default (NULL), the
  length of the corresponding query (\`seq\`) is assigned for every
  record.

- seq:

  character vector of query sequences. When default (NULL), random
  sequence is assigned. The lengths of these random sequences equal to
  the lengths of methylation call strings from the \`XM\` optional
  parameter (if supplied), or to the \`tlen\` parameter (if defined). If
  none of these parameters is supplied, length of every \`seq\` will
  equal 10.

- qual:

  query sequence quality strings (ASCII of base QUALity plus 33). When
  default (NULL), quality of every base is assigned to "F" (QUALity of
  47 + 33). The lengths of these quality strings equal to the length of
  the corresponding query sequences (\`seq\`) for every record.

- ...:

  optional tags to add to the records, in the form \`tag=value\`. Value
  can be either:

  - an integer vector to create a tag with a single integer value per
    alignment record (e.g., "NM" tag),

  - or a float vector to create a tag with a single float value per
    alignment record,

  - or a character vector (e.g., "XM" tag for methylation call string,
    "XG"/"YD"/"ZS" tag for reference strand read was aligned to)

  - or a list of numeric vectors to create tags array holding arrays of
    numeric values.

- verbose:

  boolean to report progress and timings (default: TRUE).

## Value

number of BAM records written (if \`output.bam.file\` is not NULL) or
[`data.table`](https://rdatatable.gitlab.io/data.table/reference/data.table.html)
object containing final records prepared for writing. NB: this object
has 0-based coordinates and numerically encoded reference names.

## Details

The function creates sample alignment records and saves them in BAM
file. Output can be used to test epialleleR methods as well as other
tools for methylation analysis. This method can significantly simplify
calculation of methylation metrics on example data (beta, VEF, and lMHL
values of epialleleR; methylation heterogeneity metrics of other tools).

The number of records written will be equal to the largest length of any
supplied (nondefault) parameter or 1 if no parameters were supplied. If
lengths of supplied parameters differ, shorter vectors will be recycled
(a whole number of times or with remainder if necessary).

Please note that function performs almost no validity checks for
supplied fields. In particular, be extra careful constructing paired-end
BAM alignments, and if necessary use \`samtools\` to perform validity
check or manual editing after BAM-\>SAM conversion.

## See also

[`generateCytosineReport`](generateCytosineReport.md) and
[`generateMhlReport`](generateMhlReport.md) for methylation reports at
the level of individual cytosines, as well as \`epialleleR\` vignettes
for the description of usage and sample data.

[Samtools](https://www.htslib.org/) for viewing BAM files.
[SAMv1](http://samtools.github.io/hts-specs/SAMv1.pdf) file format
specifications. Specifications of [optional SAM
tags](https://samtools.github.io/hts-specs/SAMtags.pdf).
[metheor](https://doi.org/10.1371/journal.pcbi.1010946) for ultrafast
DNA methylation heterogeneity calculation from bisulfite alignments.

## Examples

``` r
  out.bam <- tempfile(pattern="simulated", fileext=".bam")
  simulateBam(
    output.bam.file=out.bam,
    pos=c(1, 2),
    XM=c("ZZZzzZZZ", "ZZzzzzZZ"),
    XG=c("CT", "AG"),
    xi=5:6,
    xf=0.05,
    ai=list(as.integer(c(1:3)), as.integer(c(4:6))),
    af=list(seq(-1, 1, 0.5))
  )
#> Writing sample BAM 
#> [0.003s]
#> [1] 2
  generateCytosineReport(out.bam, threshold.reads=FALSE)
#> Checking BAM file: 
#> short-read, single-end, unsorted alignment detected
#> Reading single-end BAM file 
#> [0.002s]
#> Filtering reads 
#> [0.000s]
#> Preparing cytosine report 
#> [0.001s]
#>      rname strand   pos context  meth unmeth
#>     <fctr> <fctr> <int>  <fctr> <int>  <int>
#>  1:   chrS      +     1      CG     1      0
#>  2:   chrS      +     2      CG     1      0
#>  3:   chrS      -     2      CG     1      0
#>  4:   chrS      +     3      CG     1      0
#>  5:   chrS      -     3      CG     1      0
#>  6:   chrS      +     4      CG     0      1
#>  7:   chrS      -     4      CG     0      1
#>  8:   chrS      +     5      CG     0      1
#>  9:   chrS      -     5      CG     0      1
#> 10:   chrS      +     6      CG     1      0
#> 11:   chrS      -     6      CG     0      1
#> 12:   chrS      +     7      CG     1      0
#> 13:   chrS      -     7      CG     0      1
#> 14:   chrS      +     8      CG     1      0
#> 15:   chrS      -     8      CG     1      0
#> 16:   chrS      -     9      CG     1      0
  # check this BAM with `samtools view` or using `output.bam.file=NULL`
```

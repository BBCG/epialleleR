# generateBedReport

\`generateBedReport\`, \`generateAmpliconReport\`,
\`generateCaptureReport\` – these functions match BAM reads to the set
of genomic locations and return the fraction of reads with an average
methylation level passing an arbitrary threshold.

## Usage

``` r
generateAmpliconReport(
  bam,
  bed,
  report.file = NULL,
  zero.based.bed = FALSE,
  match.tolerance = 1,
  threshold.reads = TRUE,
  threshold.context = c("CG", "CHG", "CHH", "CxG", "CX"),
  min.context.sites = 2,
  min.context.beta = 0.5,
  max.outofcontext.beta = 0.1,
  ...,
  gzip = FALSE,
  verbose = TRUE
)

generateCaptureReport(
  bam,
  bed,
  report.file = NULL,
  zero.based.bed = FALSE,
  match.min.overlap = 1,
  threshold.reads = TRUE,
  threshold.context = c("CG", "CHG", "CHH", "CxG", "CX"),
  min.context.sites = 2,
  min.context.beta = 0.5,
  max.outofcontext.beta = 0.1,
  ...,
  gzip = FALSE,
  verbose = TRUE
)

generateBedReport(
  bam,
  bed,
  report.file = NULL,
  zero.based.bed = FALSE,
  bed.type = c("amplicon", "capture"),
  match.tolerance = 1,
  match.min.overlap = 1,
  threshold.reads = TRUE,
  threshold.context = c("CG", "CHG", "CHH", "CxG", "CX"),
  min.context.sites = 2,
  min.context.beta = 0.5,
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

- bed:

  Browser Extensible Data (BED) file location string OR object of class
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  holding genomic coordinates for regions of interest. The style of
  seqlevels of BED file/object must be the same as the style of
  seqlevels of BAM file/object used. The
  BED/[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  rows are **not** sorted internally. As of now, the strand information
  is ignored and reads (read pairs) matching both strands are separately
  counted and reported.

- report.file:

  file location string to write the BED report. If NULL (the default)
  then report is returned as a
  [`data.table`](https://rdatatable.gitlab.io/data.table/reference/data.table.html)
  object.

- zero.based.bed:

  boolean defining if BED coordinates are zero based (default: FALSE).

- match.tolerance:

  integer for the largest difference between read's and BED
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  start or end positions during matching of amplicon-based NGS reads
  (default: 1).

- threshold.reads:

  boolean defining if sequence reads should be thresholded before
  counting reads belonging to variant epialleles (default: TRUE).
  Disabling thresholding is possible but makes no sense in the context
  of this function, because all the reads will be assigned to the
  variant epiallele, which will result in VEF==1 (in such case \`NA\`
  VEF values are returned in order to avoid confusion). As thresholding
  is **not** recommended for long-read sequencing data, this function is
  **not** recommended for such data either.

- threshold.context:

  string defining cytosine methylation context used for thresholding the
  reads:

  - "CG" (the default) – within-the-context: CpG cytosines (called as
    zZ), out-of-context: all the other cytosines (hHxX)

  - "CHG" – within-the-context: CHG cytosines (xX), out-of-context: hHzZ

  - "CHH" – within-the-context: CHH cytosines (hH), out-of-context: xXzZ

  - "CxG" – within-the-context: CG and CHG cytosines (zZxX),
    out-of-context: CHH cytosines (hH)

  - "CX" – all cytosines are considered within-the-context, this
    effectively results in no thresholding

  This option has no effect when read thresholding is disabled.

- min.context.sites:

  non-negative integer for minimum number of cytosines within the
  \`threshold.context\` (default: 2). Reads containing **fewer**
  within-the-context cytosines are considered completely unmethylated
  (thus belonging to the reference epiallele). This option has no effect
  when read thresholding is disabled.

- min.context.beta:

  real number in the range \[0;1\] (default: 0.5). Reads with average
  beta value for within-the-context cytosines **below** this threshold
  are considered completely unmethylated (thus belonging to the
  reference epiallele). This option has no effect when read thresholding
  is disabled.

- max.outofcontext.beta:

  real number in the range \[0;1\] (default: 0.1). Reads with average
  beta value for out-of-context cytosines **above** this threshold are
  considered completely unmethylated (thus belonging to the reference
  epiallele). This option has no effect when read thresholding is
  disabled.

- ...:

  other parameters to pass to the [`preprocessBam`](preprocessBam.md)
  function. Options have no effect if preprocessed BAM data was supplied
  as an input.

- gzip:

  boolean to compress the report (default: FALSE).

- verbose:

  boolean to report progress and timings (default: TRUE).

- match.min.overlap:

  integer for the smallest overlap between read's and BED
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  start or end positions during matching of capture-based NGS reads
  (default: 1). If read matches two or more BED genomic regions, only
  the first match is taken (input
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  are **not** sorted internally).

- bed.type:

  character string for the type of assay that was used to produce
  sequencing reads:

  - "amplicon" (the default) – used for amplicon-based next-generation
    sequencing when exact coordinates of sequenced fragments are known.
    Matching of reads to genomic ranges are then performed by the read's
    start or end positions, either of which should be no further than
    \`match.tolerance\` bases away from the start or end position of
    genomic ranges given in BED
    file/[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
    object

  - "capture" – used for capture-based next-generation sequencing when
    reads partially overlap with the capture target regions. Read is
    considered to match the genomic range when their overlap is more or
    equal to \`match.min.overlap\`. If read matches two or more BED
    genomic regions, only the first match is taken (input
    [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
    are **not** sorted internally)

## Value

[`data.table`](https://rdatatable.gitlab.io/data.table/reference/data.table.html)
object containing VEF report for BED
[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html) or
NULL if report.file was specified. If BAM file contains reads that would
not match to any of BED
[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html),
the last row in the report will contain information on such reads (with
seqnames, start and end equal to NA). The report columns are:

- seqnames – reference sequence name

- start – start of genomic region

- end – end of genomic region

- width – width of genomic region

- strand – strand

- ... – other columns that were present in BED or metadata columns of
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object

- nreads+ – number of reads (pairs) mapped to the forward ("+") strand

- nreads- – number of reads (pairs) mapped to the reverse ("-") strand

- VEF – frequency of reads passing the threshold

## Details

Functions report hypermethylated variant epiallele frequencies (VEF) per
genomic region of interest using BAM and BED files as input. Reads (for
paired-end sequencing alignment files - read pairs as a single entity)
are matched to genomic locations by exact coordinates
(\`generateAmpliconReport\` or \`generateBedReport\` with an option
bed.type="amplicon") or minimum overlap (\`generateCaptureReport\` or
\`generateBedReport\` with an option bed.type="capture") – the former to
be used for amplicon-based NGS data, while the latter – for the
capture-based NGS data. The function's logic is explained below.

Let's suppose we have a BAM file with four reads, all mapped to the "+"
strand of chromosome 1, positions 1-16. The genomic range is supplied as
a parameter \`bed = as("chr1:1-100", "GRanges")\`. Assuming the default
values for the thresholding parameters (threshold.reads = TRUE,
threshold.context = "CG", min.context.sites = 2, min.context.beta = 0.5,
max.outofcontext.beta = 0.1), the input and results will look as
following:

|  |  |  |
|----|----|----|
| methylation string | threshold | explained |
| ...Z..x+.h..x..h. | below | min.context.sites \< 2 (only one zZ base) |
| ...Z..z.h..x..h. | above | pass all criteria |
| ...Z..z.h..X..h. | below | max.outofcontext.beta \> 0.1 (1XH / 3xXhH = 0.33) |
| ...Z..z.h..z-.h. | below | min.context.beta \< 0.5 (1Z / 3zZ = 0.33) |

Only the second read will satisfy all of the thresholding criteria,
leading to the following BED report (given that all reads map to
chr1:+:1-16):

|          |       |     |       |        |         |         |      |
|----------|-------|-----|-------|--------|---------|---------|------|
| seqnames | start | end | width | strand | nreads+ | nreads- | VEF  |
| chr1     | 1     | 100 | 100   | \*     | 4       | 0       | 0.25 |

Please note, that read thresholding by an average methylation level (as
explained above) makes little sense for long-read sequencing alignments,
as such reads can cover multiple regions with very different DNA
methylation properties. Instead, please use
[`extractPatterns`](extractPatterns.md), limiting pattern output to the
region of interest only.

## See also

[`preprocessBam`](preprocessBam.md) for preloading BAM data,
[`generateCytosineReport`](generateCytosineReport.md) for methylation
statistics at the level of individual cytosines,
[`generateVcfReport`](generateVcfReport.md) for evaluating epiallele-SNV
associations, [`extractPatterns`](extractPatterns.md) for exploring
methylation patterns and [`plotPatterns`](plotPatterns.md) for pretty
plotting of its output, [`generateBedEcdf`](generateBedEcdf.md) for
analysing the distribution of per-read beta values, and \`epialleleR\`
vignettes for the description of usage and sample data.

[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
class for working with genomic ranges,
[`seqlevelsStyle`](https://rdrr.io/pkg/GenomeInfoDb/man/seqlevelsStyle.html)
function for getting or setting the seqlevels style.

## Examples

``` r
  # amplicon data
  amplicon.bam    <- system.file("extdata", "amplicon010meth.bam",
                                 package="epialleleR")
  amplicon.bed    <- system.file("extdata", "amplicon.bed",
                                 package="epialleleR")
  amplicon.report <- generateAmpliconReport(bam=amplicon.bam,
                                            bed=amplicon.bed)
#> Reading BED file 
#> [0.007s]
#> Checking BAM file: 
#> short-read, paired-end, name-sorted alignment detected
#> Reading paired-end BAM file 
#> [0.004s]
#> Thresholding reads 
#> [0.000s]
#> Preparing amplicon report 
#> [0.013s]
  
  # capture NGS
  capture.bam    <- system.file("extdata", "capture.bam",
                                package="epialleleR")
  capture.bed    <- system.file("extdata", "capture.bed",
                                package="epialleleR")
  capture.report <- generateCaptureReport(bam=capture.bam, bed=capture.bed)
#> Reading BED file 
#> [0.014s]
#> Checking BAM file: 
#> short-read, paired-end, name-sorted alignment detected
#> Reading paired-end BAM file 
#> [0.013s]
#> Thresholding reads 
#> [0.002s]
#> Preparing capture report 
#> [0.018s]
  
  # generateAmpliconReport and generateCaptureReport are just aliases
  # of the generateBedReport
  bed.report <- generateBedReport(bam=capture.bam, bed=capture.bed,
                                  bed.type="capture")
#> Reading BED file 
#> [0.009s]
#> Checking BAM file: 
#> short-read, paired-end, name-sorted alignment detected
#> Reading paired-end BAM file 
#> [0.013s]
#> Thresholding reads 
#> [0.001s]
#> Preparing capture report 
#> [0.017s]
  identical(capture.report, bed.report)
#> [1] TRUE
```

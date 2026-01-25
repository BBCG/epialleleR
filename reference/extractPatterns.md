# extractPatterns

This function extracts methylation patterns (epialleles) for a given
genomic region of interest.

## Usage

``` r
extractPatterns(
  bam,
  bed,
  bed.row = 1,
  zero.based.bed = FALSE,
  match.min.overlap = 1,
  extract.context = c("CG", "CHG", "CHH", "CxG", "CX"),
  min.context.freq = 0.01,
  clip.patterns = FALSE,
  strand.offset = c(CG = 1, CHG = 2, CHH = 0, CxG = 0, CX = 0)[extract.context],
  highlight.positions = c(),
  ...,
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
  holding genomic coordinates for regions of interest. It is used to
  match sequencing reads to the genomic regions for pattern extraction.
  The style of seqlevels of BED file/object must match the style of
  seqlevels of the BAM file/object used. The
  BED/[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  rows are **not** sorted internally. As of now, the strand information
  is ignored and patterns matching both strands are extracted.

- bed.row:

  single non-negative integer specifying what \`bed\` region should be
  included in the output (default: 1).

- zero.based.bed:

  boolean defining if BED coordinates are zero based (default: FALSE).

- match.min.overlap:

  integer for the smallest overlap between read's and
  BED/[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  start or end positions during matching of capture-based NGS reads
  (default: 1).

- extract.context:

  string defining cytosine methylation context used to report:

  - "CG" (the default) – CpG cytosines (called as zZ)

  - "CHG" – CHG cytosines (xX)

  - "CHH" – CHH cytosines (hH)

  - "CxG" – CG and CHG cytosines (zZxX)

  - "CX" – all cytosines

- min.context.freq:

  real number in the range \[0;1\] (default: 0.01). Genomic positions
  that are covered by smaller fraction of patterns (e.g., with erroneous
  context) won't be included in the report.

- clip.patterns:

  boolean if patterns should not extend over the edge of \`bed\` region
  (default: FALSE).

- strand.offset:

  single non-negative integer specifying the offset of bases at the
  reverse (-) strand compared to the forward (+) strand. Allows to
  "merge" genomic positions when methylation is symmetric (in CG and CHG
  contexts). By default, equals 1 for \`extract.context\`=="CG", 2 for
  "CHG", or 0 otherwise.

- highlight.positions:

  integer vector with genomic positions of bases to include in every
  overlapping pattern. Allows to visualize the distribution of
  single-nucleotide variations (SNVs) among methylation patterns.
  \`highlight.positions\` takes precedence if any of these positions
  overlap with within-the-context positions of methylation pattern.

- ...:

  other parameters to pass to the [`preprocessBam`](preprocessBam.md)
  function. Options have no effect if preprocessed BAM data was supplied
  as an input.

- verbose:

  boolean to report progress and timings (default: TRUE).

## Value

[`data.table`](https://rdrr.io/pkg/data.table/man/data.table.html)
object containing per-read (pair) base methylation information for the
genomic region of interest. The report columns are:

- seqnames – read (pair) reference sequence name

- strand – read (pair) strand

- start – start of the read (pair)

- end – end of the read (pair)

- nbase – number of within-the-context bases for this read (pair)

- beta – beta value of this read (pair), calculated as a ratio of the
  number of methylated within-the-context bases to the total number of
  within-the-context bases

- pattern – hex representation of 64-bit FNV-1a hash calculated for all
  reported base positions and bases in this read (pair). This hash value
  depends only on included genomic positions and their methylation call
  string chars (hHxXzZ) or nucleotides (ACGT, for highlighted bases
  only), thus it is expected to be unique for every methylation pattern,
  although equal for identical methylation patterns independently on
  read (pair) start, end, or strand (when correct \`strand.offset\` is
  given)

- ... – columns for each genomic position that hold corresponding
  methylation call string char, or NA if position is not present in the
  read (pair)

## Details

The function matches reads (for paired-end sequencing alignment files -
read pairs as a single entity) to the genomic region provided in a BED
file/[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object, extracts methylation statuses of bases within those reads, and
returns a data frame which can be used for further analysis and/or
plotting of DNA methylation patterns by
[`plotPatterns`](plotPatterns.md) function.

## See also

[`plotPatterns`](plotPatterns.md) for pretty plotting of the output,
[`preprocessBam`](preprocessBam.md) for preloading BAM data,
[`generateCytosineReport`](generateCytosineReport.md) for methylation
statistics at the level of individual cytosines,
[`generateBedReport`](generateBedReport.md) for genomic region-based
statistics, [`generateVcfReport`](generateVcfReport.md) for evaluating
epiallele-SNV associations, [`generateBedEcdf`](generateBedEcdf.md) for
analysing the distribution of per-read beta values, and \`epialleleR\`
vignettes for the description of usage and sample data.

## Examples

``` r
  # amplicon data
  amplicon.bam <- system.file("extdata", "amplicon010meth.bam",
                              package="epialleleR")
  amplicon.bed <- system.file("extdata", "amplicon.bed",
                              package="epialleleR")
  
  # extract patterns
  patterns <- extractPatterns(bam=amplicon.bam, bed=amplicon.bed, bed.row=3)
#> Reading BED file 
#> [0.032s]
#> Checking BAM file: 
#> short-read, paired-end, name-sorted alignment detected
#> Reading paired-end BAM file 
#> [0.009s]
#> Extracting methylation patterns 
#> [0.026s]
  
  # and then plot them
  plotPatterns(patterns)
#> 238 patterns supplied
#> 45 unique
#> 9 most frequent unique patterns were selected for plotting using 10 beta value bins:
#> [0,0.1) [0.1,0.2) [0.2,0.3) [0.3,0.4) [0.4,0.5) [0.5,0.6) [0.6,0.7) [0.7,0.8) [0.8,0.9) [0.9,1]
#>       2         1         1         0         0         0         0         1         2       2

  
  # patterns from long reads, clipped to two narrow areas of interest
  long.bam <- system.file("extdata", "longread.bam", package="epialleleR")
  long.bed <- as(c("chr17:43125000-43125999", "chr17:43126001-43126999"),
                 "GRanges")
  long.data <- preprocessBam(
    bam=long.bam, targets=long.bed, clip.to.targets=TRUE,
    min.mapq=30, min.baseq=20, min.prob=178
  )
#> Checking BAM file: 
#> long-read, single-end, unsorted alignment detected
#> Reading single-end BAM file 
#> [0.028s]
  plotPatterns(
    extractPatterns(bam=long.data,
                    bed=as("chr17:43125000-43127000", "GRanges")),
    npatterns.per.bin=Inf
  )
#> Extracting methylation patterns 
#> [0.053s]
#> 37 patterns supplied
#> 37 unique
#> 37 most frequent unique patterns were selected for plotting using 10 beta value bins:
#> [0,0.1) [0.1,0.2) [0.2,0.3) [0.3,0.4) [0.4,0.5) [0.5,0.6) [0.6,0.7) [0.7,0.8) [0.8,0.9) [0.9,1]
#>      16         0         0         1         1         1         1         0         6      11

  
```

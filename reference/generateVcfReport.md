# generateVcfReport

This function reports base frequencies at particular genomic positions
and tests their association with the methylation status of the
sequencing reads.

## Usage

``` r
generateVcfReport(
  bam,
  vcf,
  vcf.style = NULL,
  bed = NULL,
  report.file = NULL,
  zero.based.bed = FALSE,
  cytosine.context = c("CG", "CHG", "CHH", "CxG", "CX"),
  filter.reads = TRUE,
  min.context.sites = 0,
  max.outofcontext.beta = 0.1,
  threshold.reads = TRUE,
  min.context.beta = 0.5,
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

- vcf:

  Variant Call Format (VCF) file location string OR a VCF object
  returned by
  [`readVcf`](https://rdrr.io/pkg/VariantAnnotation/man/readVcf-methods.html)
  function. If VCF object is supplied, the style of its seqlevels must
  match the style of seqlevels of the BAM file/object used.

- vcf.style:

  string for the seqlevels style of the VCF file, if different from BED
  file/object. Only has effect when \`vcf\` parameter points to the VCF
  file location and \`bed\` is not NULL. Possible values:

  - NULL (the default) – seqlevels in BED file/object and VCF file are
    the same

  - "NCBI", "UCSC", ... – valid parameters of
    [`seqlevelsStyle`](https://rdrr.io/pkg/GenomeInfoDb/man/seqlevelsStyle.html)
    function

- bed:

  Browser Extensible Data (BED) file location string OR object of class
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  holding genomic coordinates for regions of interest. It is used to
  include only the specific genomic ranges when the VCF file is loaded.
  This option has no effect when VCF object is supplied as a \`vcf\`
  parameter. The style of seqlevels of BED file/object must match the
  style of seqlevels of the BAM file/object used.

- report.file:

  file location string to write the VCF report. If NULL (the default)
  then report is returned as a
  [`data.table`](https://rdatatable.gitlab.io/data.table/reference/data.table.html)
  object.

- zero.based.bed:

  boolean defining if BED coordinates are zero based (default: FALSE).

- cytosine.context:

  string defining cytosine methylation context used for filtering and/or
  thresholding the reads:

  - "CG" (the default) – within-the-context: CpG cytosines (called as
    zZ), out-of-context: all the other cytosines (hHxX)

  - "CHG" – within-the-context: CHG cytosines (xX), out-of-context: hHzZ

  - "CHH" – within-the-context: CHH cytosines (hH), out-of-context: xXzZ

  - "CxG" – within-the-context: CG and CHG cytosines (zZxX),
    out-of-context: CHH cytosines (hH)

  - "CX" – all cytosines are considered within-the-context, this
    effectively results in no thresholding

- filter.reads:

  boolean defining if sequence reads with too few context bases or too
  high out-of-context cytosine methylation should be filtered out (e.g.,
  reads resulting from incompletely bisulfite-converted templates).
  Default: TRUE. Filtering is strongly recommended for short-read
  sequencing (bisulfite or enzymatic) because it removes reads from
  incompletely converted DNA molecules.

- min.context.sites:

  non-negative integer for minimum number of cytosines within the
  \`cytosine.context\` (default: 0, i.e., all reads will satisfy this
  criterion). When \`min.context.sites\`\>0, reads containing **fewer**
  within-the-context cytosines will not be thresholded and will be
  ignored in further computations. This option has no effect when read
  filtering is disabled.

- max.outofcontext.beta:

  real number in the range \[0;1\] (default: 0.1). Reads with average
  beta value for out-of-context cytosines **above** this threshold will
  not be thresholded and will be ignored in further computations. This
  option has no effect when read filtering is disabled.

- threshold.reads:

  boolean defining if sequence reads should be thresholded before
  counting bases in reference and variant epialleles (default: TRUE).
  Disabling thresholding is possible but makes no sense in the context
  of this function, because all the reads will be assigned to the
  variant epiallele, which will result in Fisher's Exact test p-value of
  1 (in columns \`FEp+\` and \`FEP-\`).

- min.context.beta:

  real number in the range \[0;1\] (default: 0.5). Reads with average
  beta value for within-the-context cytosines **below** this threshold
  are considered completely unmethylated (thus belonging to the
  reference epiallele). This option has no effect when read thresholding
  is disabled.

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
object containing VCF report or NULL if report.file was specified. The
report columns are:

- name – variation identifier (e.g. "rs123456789")

- seqnames – reference sequence name

- range – genomic coordinates of the variation

- REF – base at the reference allele

- ALT – base at the alternative allele

- nfiltered – number of filtered out reads

- \[M\|U\]\[+\|-\]\[Ref\|Alt\] – number of **Ref**erence or
  **Alt**ernative bases that were found at this particular position
  within **M**ethylated (above threshold) or **U**nmethylated (below
  threshold) reads that were mapped to **"+"** (forward) or **"-"**
  (reverse) DNA strand. NA values mean that it is not possible to
  determine the number of bases due to the bisulfite conversion-related
  limitations (C-\>T variants on "+" and G-\>A on "-" strands)

- SumRef – sum of all **Ref**erence base counts

- SumAlt – sum of all **Alt**ernative base counts

- FEp+ – Fisher Exact test p-value for association of a variation with
  methylation status of the reads that map to the **"+"** (forward) DNA
  strand. Calculated using following contingency table:

  |       |       |
  |-------|-------|
  | M+Ref | M+Alt |
  | U+Ref | U+Alt |

- FEp- – Fisher Exact test p-value for association of a variation with
  methylation status of the reads that map to the **"-"** (reverse) DNA
  strand. Calculated using following contingency table:

  |       |       |
  |-------|-------|
  | M-Ref | M-Alt |
  | U-Ref | U-Alt |

## Details

Using BAM reads and sequence variation information as an input,
\`generateVcfReport\` function filters and thresholds the reads (for
paired-end sequencing alignment files - read pairs as a single entity)
according to supplied parameters and calculates the occurrence of
**Ref**erence and **Alt**ernative bases within reads, taking into the
account DNA strand the read mapped to and average methylation level
(epiallele status) of the read.

The information on sequence variation can be supplied as a Variant Call
Format (VCF) file location or an object of class VCF, returned by the
[`readVcf`](https://rdrr.io/pkg/VariantAnnotation/man/readVcf-methods.html)
function call. As whole-genome VCF files can be extremely large, it is
strongly advised to use only relevant subset of their data, prefiltering
the VCF object manually before calling \`generateVcfReport\` or
specifying \`bed\` parameter when \`vcf\` points to the location of such
large VCF file. Please note that all the BAM, BED and VCF files must use
the same style for seqlevels (i.e. chromosome names).

After counting, function checks if certain bases occur more often within
reads belonging to certain epialleles using Fisher Exact test (HTSlib's
own implementation) and reports separate p-values for reads mapped to
**"+"** (forward) and **"-"** (reverse) DNA strands.

Please note that the final report currently includes only the VCF
entries with single-base REF and ALT alleles. Also, the default
(\`min.baseq=0\`) output of \`generateVcfReport\` is equivalent to the
one of \`samtools mplieup -Q 0 ...\`, and therefore may result in false
SNVs caused by misalignments. Remember to increase \`min.baseq\`
(\`samtools mpileup -Q\` default value is 13) to obtain higher-quality
results.

Read thresholding by an average methylation level used in this function
makes little sense for long-read sequencing alignments, as such reads
can cover multiple regions with very different DNA methylation
properties. If necessary, one could either clip long sequencing reads to
narrow \`targets\` in [`preprocessBam`](preprocessBam.md) function
during BAM loading or use [`extractPatterns`](extractPatterns.md),
limiting pattern output to the region of interest only.

## See also

[`preprocessBam`](preprocessBam.md) for preloading BAM data,
[`generateCytosineReport`](generateCytosineReport.md) for methylation
statistics at the level of individual cytosines,
[`generateBedReport`](generateBedReport.md) for genomic region-based
statistics, [`extractPatterns`](extractPatterns.md) for exploring
methylation patterns and [`plotPatterns`](plotPatterns.md) for pretty
plotting of its output, [`generateBedEcdf`](generateBedEcdf.md) for
analysing the distribution of per-read beta values, and \`epialleleR\`
vignettes for the description of usage and sample data.

[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
class for working with genomic ranges,
[`readVcf`](https://rdrr.io/pkg/VariantAnnotation/man/readVcf-methods.html)
function for loading VCF data,
[`seqlevelsStyle`](https://rdrr.io/pkg/GenomeInfoDb/man/seqlevelsStyle.html)
function for getting or setting the seqlevels style.

## Examples

``` r
  capture.bam <- system.file("extdata", "capture.bam", package="epialleleR")
  capture.bed <- system.file("extdata", "capture.bed", package="epialleleR")
  capture.vcf <- system.file("extdata", "capture.vcf.gz",
                             package="epialleleR")
  
  # VCF report
  vcf.report <- generateVcfReport(bam=capture.bam, bed=capture.bed,
                                  vcf=capture.vcf)
#> Loading required namespace: VariantAnnotation
#> Reading BED file 
#> [0.023s]
#> Reading VCF file 
#> [5.768s]
#> Checking BAM file: 
#> short-read, paired-end, name-sorted alignment detected
#> Reading paired-end BAM file 
#> [0.012s]
#> Filtering and thresholding reads 
#> [0.001s]
#> Extracting base frequences 
#> [0.146s]
  
  # toy example to illustrate the logic of computations
  if (requireNamespace("VariantAnnotation", quietly=TRUE)) {
    # simulate toy BAM
    temp.bam <- tempfile(fileext=".bam")
    simulateBam(output.bam.file=temp.bam, rname="chr1", XG="CT",
                seq=c("AGACGTTAGTAATAGTA", "AGACGTTGTAATAGTA",
                      "AAACGTTGTAACAGTA",  "AAACGTTGTAATGTA"),
                XM=c( "...Z..x+.h..x..h.", "...Z..z.h..x..h.",
                      "...Z..z.h..X..h.",  "...Z..z.h..z.h."),
                cigar=c("7M1I9M", "16M", "16M", "12M1D3M"))
    # toy VCF
    vcf <- VariantAnnotation::VCF(rowRanges=as("chr1:2", "GRanges"),
                                  collapsed=FALSE)
    VariantAnnotation::ref(vcf) <- as("A", "DNAStringSet")
    VariantAnnotation::alt(vcf) <- as("G", "DNAStringSet")
    
    # when default values of filtering and thresholding parameters are used,
    # read filtering will exclude the third read from this BAM file
    # because it has too many outside-of-context methylated cytosines.
    
    # results with read filtering and thresholding
    generateVcfReport(bam=temp.bam, vcf=vcf)
    # results without read filtering
    generateVcfReport(bam=temp.bam, vcf=vcf, filter.reads=FALSE)
  }
#> Writing sample BAM 
#> [0.003s]
#> Checking BAM file: 
#> short-read, single-end, unsorted alignment detected
#> Reading single-end BAM file 
#> [0.002s]
#> Filtering and thresholding reads 
#> [0.000s]
#> Extracting base frequences 
#> [0.033s]
#> Checking BAM file: 
#> short-read, single-end, unsorted alignment detected
#> Reading single-end BAM file 
#> [0.002s]
#> Thresholding reads 
#> [0.000s]
#> Extracting base frequences 
#> [0.033s]
#>    seqnames range    REF    ALT nfiltered M+Ref U+Ref M-Ref U-Ref M+Alt U+Alt
#>      <fctr> <int> <char> <char>    <lgcl> <num> <num> <num> <num> <num> <num>
#> 1:     chr1     2      A      G        NA     1     1    NA    NA     2     0
#>    M-Alt U-Alt SumRef SumAlt  FEp+  FEp-
#>    <num> <num>  <num>  <num> <num> <num>
#> 1:    NA    NA      2      2     1    NA
```

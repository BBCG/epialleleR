# Fast, epiallele-aware methylation caller and reporter [![logo](articles/epialleleR_logo.svg)](https://github.com/BBCG/epialleleR)

[![](https://github.com/BBCG/epialleleR/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/BBCG/epialleleR/actions)
[![](https://codecov.io/gh/BBCG/epialleleR/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/BBCG/epialleleR/tree/devel)
[![](https://bioconductor.org/shields/years-in-bioc/epialleleR.svg)](https://bioconductor.org/packages/release/bioc/html/epialleleR.html)
[![install from
r-universe](https://bioc.r-universe.dev/epialleleR/badges/version)](https://bioc.r-universe.dev/epialleleR)

## Introduction

*`epialleleR`* is an R package for calling and reporting cytosine DNA
methylation in many useful ways. Developed to help identify and quantify
epimutations (aberrant DNA methylation), it has now acquired multiple
additional functions to dissect DNA methylation. But the main feature of
the package is to report frequencies of epimutations (variant epiallele
frequencies, VEF) at the level of genomic regions or individual
cytosines. All you need in order to use it is a binary alignment map
(BAM) file from basically any next-generation (methylation or native)
sequencing experiment.

![](./articles/epialleles.png)

### Input Data

- short-read and long-read (native)
- paired-end and single-end
- whole-genome, genome-wide (e.g., hybridization capture or adaptive
  sampling), and narrowly targeted (e.g., amplicon panels)

### Current Features

- calling cytosine methylation and saving calls in a new BAM file
  (*`callMethylation`*)
- creating sample BAM files from scratch given mandatory and optional
  BAM fields (*`simulateBam`*)
- conventional reporting of cytosine methylation
  (*`generateCytosineReport`*)
- reporting the hypermethylated variant epiallele frequency (VEF) at the
  level of genomic regions (*`generate[Bed|Amplicon|Capture]Report`*) or
  individual cytosines (*`generateCytosineReport`*)
- reporting linearised Methylated Haplotype Load (lMHL,
  *`generateMhlReport`*)
- extracting methylation patterns for genomic region of interest
  (*`extractPatterns`*)
- visualising methylation patterns (*`plotPatterns`*)
- testing for the association between epiallele methylation status and
  sequence variations (*`generateVcfReport`*)
- assessing the distribution of per-read beta values for genomic regions
  of interest (*`generateBedEcdf`*)

### Recent improvements

##### v1.20 \[BioC 3.23\]

- using genomic coordinates of targets, only a subset of BAM reads or
  only fragments of BAM reads that are overlapping the targets can now
  be loaded
- disrupting changes in all `generate*Report` functions (from version
  1.19.1 onwards):
  - `cytosine.context` parameter instead of
    `threshold.context`/`haplotype.context`
  - new parameter `filter.reads` regulates filtering of reads with too
    few cytosines or presumable incomplete conversion of cytosines
  - reads with too few within-the-context cytosines (less than
    `min.context.sites`) or out-of-context cytosine methylation higher
    than `max.outofcontext.beta` are filtered out (discarded) instead of
    being counted as hypomethylated reads (as was in v1.19.0 and
    earlier)
  - new default value of 0 for `min.context.sites`

##### v1.14 \[BioC 3.20\]

- creates pretty plots of methylation patterns

##### v1.12 \[BioC 3.19\]

- inputs long-read sequencing alignments
- full support for short-read sequencing alignments by Illumina DRAGEN,
  Bismark, bwa-meth, BSMAP
- RRBS-specific options
- lower memory usage

##### v1.10 \[BioC 3.18\]

- inputs both single-end and paired-end sequencing alignments
- makes and stores methylation calls
- creates sample BAM files
- reports linearised MHL

##### v1.4 \[BioC 3.15\]

- significant speed-up
- method to extract and visualize methylation patterns

##### v1.2 \[BioC 3.14\]

- even faster and more memory-efficient BAM loading (by means of HTSlib)
- min.baseq parameter to reduce the effect of low quality bases on
  methylation or SNV calling (in v1.0 the output of
  *`generateVcfReport`* was equivalent to the one of
  `samtools mpileup -Q 0 ...`)

check out NEWS for more!

------------------------------------------------------------------------

## Installation

### install via Bioconductor

``` r

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("epialleleR")
```

### Install the latest version via install_github

``` r

library(devtools)
install_github("BBCG/epialleleR", build_vignettes=FALSE,
  repos=BiocManager::repositories(),
  dependencies=TRUE, type="source")
```

------------------------------------------------------------------------

## Using the package

Please read *`epialleleR`* vignette at [GitHub
pages](https://bbcg.github.io/epialleleR/articles/epialleleR.html) or
within the R environment:
[`vignette("epialleleR", package="epialleleR")`](articles/epialleleR.md),
or consult the function’s help pages for the extensive information on
usage, parameters and output values.

Comparison of beta, VEF and lMHL values for various use cases is given
by the [values](https://bbcg.github.io/epialleleR/articles/values.html)
vignette
([`vignette("values", package="epialleleR")`](articles/values.md))

Very brief synopsis:

``` r

library(epialleleR)

# make methylation calls if necessary
callMethylation(
  input.bam.file=system.file("extdata", "test", "dragen-se-unsort-xg.bam", package="epialleleR"),
  output.bam.file=tempfile(pattern="output-", fileext=".bam"),
  genome=system.file("extdata", "test", "reference.fasta.gz", package="epialleleR")
)

# make a sample BAM file from scratch
simulateBam(output.bam.file=tempfile(pattern="simulated-", fileext=".bam"),
            pos=c(1, 2), XM=c("ZZZzzZZZ", "ZZzzzzZZ"), XG=c("CT", "AG"))

# or use external files
amplicon.bam <- system.file("extdata", "amplicon010meth.bam",
                            package="epialleleR")
amplicon.bed <- system.file("extdata", "amplicon.bed", package="epialleleR")
amplicon.vcf <- system.file("extdata", "amplicon.vcf.gz", package="epialleleR")

# preload the data
bam.data <- preprocessBam(amplicon.bam)

# methylation patterns and their plot
patterns <- extractPatterns(bam=amplicon.bam, bed=amplicon.bed, bed.row=3)
plotPatterns(patterns)

# conventional cytosine report
cx.report <- generateCytosineReport(bam.data, filter.reads=FALSE,
                                    threshold.reads=FALSE, report.context="CX")

# CpG VEF report for individual bases
cg.vef.report <- generateCytosineReport(bam.data)

# BED-guided VEF report for genomic ranges
bed.report <- generateBedReport(bam=amplicon.bam, bed=amplicon.bed,
                                bed.type="capture")

# VCF report
vcf.report <- generateVcfReport(bam=amplicon.bam, bed=amplicon.bed,
                                vcf=amplicon.vcf, vcf.style="NCBI")

# lMHL report
mhl.report <- generateMhlReport(bam=amplicon.bam)
```

------------------------------------------------------------------------

### Citing the *`epialleleR`* package

Oleksii Nikolaienko, Per Eystein Lønning, Stian Knappskog, *epialleleR*:
an R/Bioconductor package for sensitive allele-specific methylation
analysis in NGS data. *GigaScience*, Volume 12, 2023, giad087,
<https://doi.org/10.1093/gigascience/giad087>. Data:
[GSE201690](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201690)

### Our experimental studies that use the package

Per Eystein Lonning, Oleksii Nikolaienko, Kathy Pan, Allison W. Kurian,
Hans Petter Petter Eikesdal, Mary Pettinger, Garnet L Anderson, Ross L
Prentice, Rowan T. Chlebowski, and Stian Knappskog. Constitutional
*BRCA1* methylation and risk of incident triple-negative breast cancer
and high-grade serous ovarian cancer. *JAMA Oncology* 2022.
<https://doi.org/10.1001/jamaoncol.2022.3846>

Oleksii Nikolaienko, Hans P. Eikesdal, Elisabet Ognedal, Bjørnar Gilje,
Steinar Lundgren, Egil S. Blix, Helge Espelid, Jürgen Geisler, Stephanie
Geisler, Emiel A.M. Janssen, Synnøve Yndestad, Laura Minsaas, Beryl
Leirvaag, Reidun Lillestøl, Stian Knappskog, Per E. Lønning. Prenatal
*BRCA1* epimutations contribute significantly to triple-negative breast
cancer development. *Genome Medicine* 2023.
<https://doi.org/10.1186/s13073-023-01262-8>. Data:
[GSE243966](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243966)

Oleksii Nikolaienko, Garnet L Anderson, Rowan T Chlebowski, Su Yon Jung,
Holly R Harris, Stian Knappskog, and Per E Lønning. *MGMT* epimutations
and risk of incident cancer of the colon, glioblastoma multiforme, and
diffuse large B-cell lymphomas. *Clinical Epigenetics* 2025.
<https://doi.org/10.1186/s13148-025-01835-x>

### *`epialleleR`* at Bioconductor

[release](https://bioconductor.org/packages/release/bioc/html/epialleleR.html),
[development
version](https://bioconductor.org/packages/devel/bioc/html/epialleleR.html)

------------------------------------------------------------------------

## License

Artistic License 2.0

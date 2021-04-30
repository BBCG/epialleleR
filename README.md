# Fast, epiallele-aware methylation reporter  <img align="right" src="inst/epialleleR_logo.svg">

[![](https://github.com/BBCG/epialleleR/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/BBCG/epialleleR/actions)
[![](https://codecov.io/gh/BBCG/epialleleR/branch/master/graph/badge.svg)](https://codecov.io/gh/BBCG/epialleleR)


## Introduction

DISCLAIMER: This is a work in progress, however the package is already usable,
and the obtained experimental results will soon be published. Main methods
(*`preprocessBam`*, *`generateCytosineReport`*, *`generateBedReport`*) won't
change. The *`generateVcfReport`* method will at some point be improved to
include variable-length sequence variations, while *`generateBedEcdf`* should
be considered somewhat experimental and may undergo significant changes or be
substituted with some other method in the future.

*`epialleleR`* is an R package for calling hypermethylated variant epiallele
frequencies (VEF) at the level of genomic regions or individual cytosines
in next-generation sequencing data using binary alignment map (BAM) files as
an input. Other functionality includes computing the empirical cumulative
distribution function for per-read beta values, and testing the significance
of the association between epiallele methylation status and base frequencies
at particular genomic positions (SNPs).

### Current Features

 * conventional reporting of cytosine methylation (*`generateCytosineReport`*)
 * calling the hypermethylated variant epiallele frequency (VEF) at the
 level of genomic regions (*`generateBedReport`*) or individual cytosines
 (*`generateCytosineReport`*)
 * assessing the distribution of per-read beta values for genomic regions of
 interest (*`generateBedEcdf`*)
 * testing for the association between epiallele methylation
 status and sequence variations (*`generateVcfReport`*)

### Coming soon

 * even faster and more memory-efficient BAM loading
 * min-baseq parameter to reduce the effect of low quality bases on SNV calling
 (at the moment the output of *`generateVcfReport`* is equivalent to the one of
 `samtools mpileup -Q 0 ...`)
 
-------

## Installation

### install via Bioconductor (NOT YET AVAILABLE)
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("epialleleR")
```

### Install the latest version via install_github
```r
library(devtools)
install_github("BBCG/epialleleR", build_vignettes=FALSE,
  repos=BiocManager::repositories(),
  dependencies=TRUE, type="source")
```

-------

## Using the package

Please read *`epialleleR`* vignettes
at [GitHub pages](https://bbcg.github.io/epialleleR/articles/epialleleR.html)
or within the R environment: `vignette("epialleleR", package="epialleleR")`, or
consult the function's help pages for the extensive information on usage,
parameters and output values.

Very brief synopsis:

```r
library(epialleleR)

# external files
amplicon.bam <- system.file("extdata", "amplicon010meth.bam",
                            package="epialleleR")
amplicon.bed <- system.file("extdata", "amplicon.bed", package="epialleleR")
amplicon.vcf <- system.file("extdata", "amplicon.vcf.gz", package="epialleleR")

# preloading the data
bam.data <- preprocessBam(amplicon.bam)

# CpG VEF report for individual bases
cg.vef.report <- generateCytosineReport(bam.data)

# BED-guided VEF report for genomic ranges
bed.report <- generateBedReport(bam=amplicon.bam, bed=amplicon.bed,
                                bed.type="capture")

# VCF report
vcf.report <- generateVcfReport(bam=amplicon.bam, bed=amplicon.bed,
                                vcf=amplicon.vcf, vcf.style="NCBI")
```

-------

### Citing the *`epialleleR`* package
The experimental data analysed using the package has not been published yet.
The citation information will be updated in the nearest future.

### *`epialleleR`* at Bioconductor
On its way now.

-------

License
---------
Artistic License/GPL
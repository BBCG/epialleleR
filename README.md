epialleleR
==========

# Introduction

DISCLAIMER: This is a work in progress, however the package is already usable,
and the obtained eperimantal results will soon be published. Main methods
(*`preprocessBam`*, *`generateCytosineReport`*, *`generateBedReport`*) won't
change. The *`generateVcfReport`* method will at some point be improved to
include variable-length sequence variations, while *`generateBedEcdf`* should
be considered somewhat experimental and may undergo significant changes or be
substituted with some other method in the future.

*`epialleleR`* is an R package for calling hypermethylated variant epiallele
(VEF) frequencies at the level of genomic regions or individual cytosines
in next-generation sequencing data using binary alignment map (BAM) files as
an input. Other functionality includes computing the empirical cumulative
distribution function for per-read beta values, and testing the significance
of the association between epiallele methylation status and base frequencies
at particular genomic positions (SNPs).

## Current Features

 * conventional reporting of cytosine methylation (*`generateCytosineReport`*)
 * calling the hypermethylated variant epiallele frequency (VEF) at the
 level of genomic regions (*`generateBedReport`*) or individual cytosines
 (*`generateCytosineReport`*)
 * assessing the distribution of per-read beta values for genomic regions of
 interest (*`generateBedEcdf`*)
 * testing for the association between epiallele methylation
 status and sequence variations (*`generateVcfReport`*)


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

### Citing the *`epialleleR`* package
[NOT YET AVAILABLE](https://doi.org/NOT.YET.AVAILABLE)

### *`epialleleR`* at Bioconductor
[NOT YET AVAILABLE](https://bioconductor.org/packages/devel/bioc/html/epialleleR.html)

-------

# How to Use

TBA


License
---------
Artistic License/GPL
epialleleR
==========

# Introduction

THIS IS A WORK IN PROGRESS, HOWEVER THE PACKAGE IS ALREADY FAIRLY USABLE. MAIN
METHODS (*`preprocessBam`*, *`generateCytosineReport`*, *`generateBedReport`*)
WON'T CHANGE. THE *`generateVcfReport`* METHOD WILL BE IMPROVED, WHILE
*`generateBedEcdf`* SHOULD BE CONSIDERED EXPERIMENTAL AND MAY DISAPPEAR OR
UNDERGO SIGNIFICANT CHANGES IN THE FUTURE.

*`epialleleR`* is an R package for calling hypermethylated
epiallele frequencies at the level of genomic regions or individual cytosines
in next-generation sequencing data using binary alignment map (BAM) files as
an input. Other functionality includes computing the empirical cumulative
distribution function for per-read beta values, and testing the significance
of the association between epiallele methylation status and base frequencies
at particular genomic positions (SNPs).

## Current Features

 * generation of cytosine report
 * generation of hypermethylated variant epiallele frequency report at the
 level of genomic regions or individual cytosines
 * computing the empirical cumulative distribution function for per-read beta
 values
 * testing the significance of the association between epiallele methylation
 status and sequence variations


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
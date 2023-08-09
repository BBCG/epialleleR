# Fast, epiallele-aware methylation <img style="float: right;" width="50px" src="inst/epialleleR_logo.svg"> caller and reporter

[![](https://github.com/BBCG/epialleleR/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/BBCG/epialleleR/actions)
[![](https://codecov.io/gh/BBCG/epialleleR/branch/devel/graph/badge.svg)](https://codecov.io/gh/BBCG/epialleleR)
[![](https://bioconductor.org/shields/years-in-bioc/epialleleR.svg)](https://bioconductor.org/packages/release/bioc/html/epialleleR.html)

## Introduction

DISCLAIMER: This is a work in progress, however the package is already usable,
and the obtained experimental results have been published. Main methods
(*`preprocessBam`*, *`generateCytosineReport`*, *`generateBedReport`*) won't
change. The *`generateVcfReport`* method might at some point be improved to
include variable-length sequence variations, while *`generateBedEcdf`* should
be considered somewhat experimental and may undergo significant changes or be
substituted with some other method in the future.

*`epialleleR`* is an R package for calling hypermethylated variant epiallele
frequencies (VEF) at the level of genomic regions or individual cytosines
in next-generation sequencing data using binary alignment map (BAM) files as
an input. Other functionality includes extracting methylation patterns,
computing the empirical cumulative distribution function for per-read beta
values, and testing the significance of the association between epiallele
methylation status and base frequencies at particular genomic positions (SNPs).

### Current Features

 * calling cytosine methylation and saving calls in BAM file
 (*`callMethylation`*)
 * conventional reporting of cytosine methylation (*`generateCytosineReport`*)
 * calling the hypermethylated variant epiallele frequency (VEF) at the
 level of genomic regions (*`generate[Bed|Amplicon|Capture]Report`*) or
 individual cytosines (*`generateCytosineReport`*)
 * extracting methylation patterns for genomic region of interest
 (*`extractPatterns`*)
 * assessing the distribution of per-read beta values for genomic regions of
 interest (*`generateBedEcdf`*)
 * testing for the association between epiallele methylation
 status and sequence variations (*`generateVcfReport`*)

### Recent improvements

##### v1.10 [BioC 3.18]

 * inputs both single-end and paired-end sequencing alignments
 * makes and stores methylation calls

 
##### v1.4 [BioC 3.15]

 * significant speed-up
 * method to visualize methylation patterns

##### v1.2 [BioC 3.14]

 * even faster and more memory-efficient BAM loading (by means of HTSlib)
 * min.baseq parameter to reduce the effect of low quality bases on 
 methylation or SNV calling (in v1.0 the output of *`generateVcfReport`* was
 equivalent to the one of `samtools mpileup -Q 0 ...`)

check out NEWS for more!
 
-------

## Installation

### install via Bioconductor
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

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

# make methylation calls if necessary
callMethylation(
  input.bam.file=system.file("extdata", "test", "dragen-se-unsort-xg.bam", package="epialleleR"),
  output.bam.file=tempfile(pattern="output-", fileext=".bam"),
  genome=system.file("extdata", "test", "reference.fasta.gz", package="epialleleR")
)

# external files
amplicon.bam <- system.file("extdata", "amplicon010meth.bam",
                            package="epialleleR")
amplicon.bed <- system.file("extdata", "amplicon.bed", package="epialleleR")
amplicon.vcf <- system.file("extdata", "amplicon.vcf.gz", package="epialleleR")

# preloading the data
bam.data <- preprocessBam(amplicon.bam)

# methylation patterns, check vignettes or method description for plotting them
patterns <- extractPatterns(bam=amplicon.bam, bed=amplicon.bed, bed.row=3)

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
Oleksii Nikolaienko, Per Eystein Lønning, Stian Knappskog, *epialleleR*: an R/Bioconductor package for sensitive allele-specific methylation analysis in NGS data. *bioRxiv* 2022.06.30.498213. [https://www.biorxiv.org/content/10.1101/2022.06.30.498213](https://www.biorxiv.org/content/10.1101/2022.06.30.498213)

### The experimental data analysed using the package
Per Eystein Lonning, Oleksii Nikolaienko, Kathy Pan, Allison W. Kurian, Hans Petter Petter Eikesdal, Mary Pettinger, Garnet L Anderson, Ross L Prentice, Rowan T. Chlebowski, and Stian Knappskog. Constitutional *BRCA1* methylation and risk of incident triple-negative breast cancer and high-grade serous ovarian cancer. *JAMA Oncology* 2022. [https://doi.org/10.1001/jamaoncol.2022.3846](https://doi.org/10.1001/jamaoncol.2022.3846)

Oleksii Nikolaienko, Hans P. Eikesdal, Bjørnar Gilje, Steinar Lundgren, Egil S. Blix, Helge Espelid, Jürgen Geisler, Stephanie Geisler, Emiel A.M. Janssen, Synnøve Yndestad, Laura Minsaas, Beryl Leirvaag, Reidun Lillestøl, Stian Knappskog, Per E. Lønning. Prenatal *BRCA1* epimutations contribute significantly to triple-negative breast cancer development. *medRxiv* 2023.05.14.23289949. [https://www.medrxiv.org/content/10.1101/2023.05.14.23289949](https://www.medrxiv.org/content/10.1101/2023.05.14.23289949)

### *`epialleleR`* at Bioconductor
[release](https://bioconductor.org/packages/release/bioc/html/epialleleR.html), 
[development version](https://bioconductor.org/packages/devel/bioc/html/epialleleR.html)

-------

License
---------
Artistic License 2.0

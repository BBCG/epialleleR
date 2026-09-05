# The epialleleR output values

Abstract

A comparison and visualisation of epialleleR output values for various
input files

## Introduction

The best possible explanation on VEF and lMHL values is given in help
files for *`generateCytosineReport`* and *`generateMhlReport`* methods,
respectively. Here we try to show some simplified and real situations,
i.e., different methylation patterns that may exist, and provide a
visual summary of *`epialleleR`* output.

The readers are welcome to try their own real and simulated data. If it
might be of interest to others, please create an issue and these
examples might get included in this vignette.

NB: the `plotMetrics` function used below is a piece of spaghetti code,
hence hidden. If you still want to use it or see what it does - browse a
[source
code](https://github.com/BBCG/epialleleR/blob/devel/vignettes/values.Rmd)
of this vignette online.

`out.bam`` ``<-`` `[`tempfile`](https://rdrr.io/r/base/tempfile.html)`(``pattern``=``"simulated"``, fileext``=``".bam"``)`` `[`set.seed`](https://rdrr.io/r/base/Random.html)`(``1``)`` `` ``# no epimutations`` `[`simulateBam`](../reference/simulateBam.md)`(`` `` output.bam.file``=``out.bam``,`` `` XM``=`[`c`](https://rdrr.io/r/base/c.html)`(`` `` `[`sapply`](https://rdrr.io/r/base/lapply.html)`(`` `` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(``1``:``1000``, ``function`` ``(``x``)`` `[`sample`](https://rdrr.io/r/base/sample.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``"Z"``,`[`rep`](https://rdrr.io/r/base/rep.html)`(``"z"``, ``9``)``)``, ``10``)``)``,`` `` ``paste``, collapse``=``""`` `` ``)`` `` ``)``,`` `` XG``=``"CT"`` ``)`` ``#> Writing sample BAM [0.157s]`` ``#> [1] 1000`` ``plotMetrics``(``out.bam``, ``as``(``"chrS:1-10"``, ``"GRanges"``)``, ``0``, title``=``"no epimutations"``)`

![](values_files/figure-html/unnamed-chunk-3-1.png)

` ``# one complete epimutation`` `[`simulateBam`](../reference/simulateBam.md)`(`` `` output.bam.file``=``out.bam``,`` `` XM``=`[`c`](https://rdrr.io/r/base/c.html)`(`` `` `[`paste`](https://rdrr.io/r/base/paste.html)`(`[`rep`](https://rdrr.io/r/base/rep.html)`(``"Z"``, ``10``)``, collapse``=``""``)``,`` `` `[`sapply`](https://rdrr.io/r/base/lapply.html)`(`` `` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(``1``:``999``, ``function`` ``(``x``)`` `[`sample`](https://rdrr.io/r/base/sample.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``"Z"``,`[`rep`](https://rdrr.io/r/base/rep.html)`(``"z"``, ``9``)``)``, ``10``)``)``,`` `` ``paste``, collapse``=``""`` `` ``)`` `` ``)``,`` `` XG``=``"CT"`` ``)`` ``#> Writing sample BAM [0.111s]`` ``#> [1] 1000`` ``plotMetrics``(``out.bam``, ``as``(``"chrS:1-10"``, ``"GRanges"``)``, title``=``"one complete epimutation"``)`

![](values_files/figure-html/unnamed-chunk-3-2.png)

` ``# one partial epimutation`` `[`simulateBam`](../reference/simulateBam.md)`(`` `` output.bam.file``=``out.bam``,`` `` XM``=`[`c`](https://rdrr.io/r/base/c.html)`(`` `` `[`paste`](https://rdrr.io/r/base/paste.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(`[`rep`](https://rdrr.io/r/base/rep.html)`(``"Z"``, ``4``)``, ``"z"``, ``"z"``, `[`rep`](https://rdrr.io/r/base/rep.html)`(``"Z"``, ``4``)``)``, collapse``=``""``)``,`` `` `[`sapply`](https://rdrr.io/r/base/lapply.html)`(`` `` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(``1``:``999``, ``function`` ``(``x``)`` `[`sample`](https://rdrr.io/r/base/sample.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``"Z"``,`[`rep`](https://rdrr.io/r/base/rep.html)`(``"z"``, ``9``)``)``, ``10``)``)``,`` `` ``paste``, collapse``=``""`` `` ``)`` `` ``)``,`` `` XG``=``"CT"`` ``)`` ``#> Writing sample BAM [0.110s]`` ``#> [1] 1000`` ``plotMetrics``(``out.bam``, ``as``(``"chrS:1-10"``, ``"GRanges"``)``, title``=``"one partial epimutation"``)`

![](values_files/figure-html/unnamed-chunk-3-3.png)

` ``# another partial epimutation`` `[`simulateBam`](../reference/simulateBam.md)`(`` `` output.bam.file``=``out.bam``,`` `` XM``=`[`c`](https://rdrr.io/r/base/c.html)`(`` `` ``"zZZZZZZZzz"``,`` `` `[`sapply`](https://rdrr.io/r/base/lapply.html)`(`` `` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(``1``:``999``, ``function`` ``(``x``)`` `[`sample`](https://rdrr.io/r/base/sample.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``"Z"``,`[`rep`](https://rdrr.io/r/base/rep.html)`(``"z"``, ``9``)``)``, ``10``)``)``,`` `` ``paste``, collapse``=``""`` `` ``)`` `` ``)``,`` `` XG``=``"CT"`` ``)`` ``#> Writing sample BAM [0.102s]`` ``#> [1] 1000`` ``plotMetrics``(``out.bam``, ``as``(``"chrS:1-10"``, ``"GRanges"``)``, title``=``"another partial epimutation"``)`

![](values_files/figure-html/unnamed-chunk-3-4.png)

` ``# several partial epimutations`` `[`simulateBam`](../reference/simulateBam.md)`(`` `` output.bam.file``=``out.bam``,`` `` XM``=`[`c`](https://rdrr.io/r/base/c.html)`(`` `` `[`sapply`](https://rdrr.io/r/base/lapply.html)`(`` `` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(``1``:``10``, ``function`` ``(``x``)`` `[`c`](https://rdrr.io/r/base/c.html)`(`[`rep`](https://rdrr.io/r/base/rep.html)`(``"Z"``, ``6``)``, `[`rep`](https://rdrr.io/r/base/rep.html)`(``"z"``, ``4``)``)``)``,`` `` ``paste``, collapse``=``""`` `` ``)``,`` `` `[`sapply`](https://rdrr.io/r/base/lapply.html)`(`` `` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(``1``:``999``, ``function`` ``(``x``)`` `[`sample`](https://rdrr.io/r/base/sample.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``"Z"``,`[`rep`](https://rdrr.io/r/base/rep.html)`(``"z"``, ``9``)``)``, ``10``)``)``,`` `` ``paste``, collapse``=``""`` `` ``)`` `` ``)``,`` `` XG``=``"CT"`` ``)`` ``#> Writing sample BAM [0.108s]`` ``#> [1] 1009`` ``plotMetrics``(``out.bam``, ``as``(``"chrS:1-10"``, ``"GRanges"``)``, title``=``"several partial epimutations"``)`

![](values_files/figure-html/unnamed-chunk-3-5.png)

` ``# several short partial epimutations`` `[`simulateBam`](../reference/simulateBam.md)`(`` `` output.bam.file``=``out.bam``,`` `` XM``=`[`c`](https://rdrr.io/r/base/c.html)`(`` `` `[`sapply`](https://rdrr.io/r/base/lapply.html)`(`` `` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(``1``:``10``, ``function`` ``(``x``)`` `[`c`](https://rdrr.io/r/base/c.html)`(`[`rep`](https://rdrr.io/r/base/rep.html)`(``"Z"``, ``4``)``, `[`rep`](https://rdrr.io/r/base/rep.html)`(``"z"``, ``6``)``)``)``,`` `` ``paste``, collapse``=``""`` `` ``)``,`` `` `[`sapply`](https://rdrr.io/r/base/lapply.html)`(`` `` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(``1``:``999``, ``function`` ``(``x``)`` `[`sample`](https://rdrr.io/r/base/sample.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``"Z"``,`[`rep`](https://rdrr.io/r/base/rep.html)`(``"z"``, ``9``)``)``, ``10``)``)``,`` `` ``paste``, collapse``=``""`` `` ``)`` `` ``)``,`` `` XG``=``"CT"`` ``)`` ``#> Writing sample BAM [0.120s]`` ``#> [1] 1009`` ``plotMetrics``(``out.bam``, ``as``(``"chrS:1-10"``, ``"GRanges"``)``, title``=``"several short partial epimutations"``)`

![](values_files/figure-html/unnamed-chunk-3-6.png)

` ``# several overlapping partial epimutations`` `[`simulateBam`](../reference/simulateBam.md)`(`` `` output.bam.file``=``out.bam``,`` `` pos``=``1``:``10``,`` `` XM``=`[`c`](https://rdrr.io/r/base/c.html)`(`` `` ``"ZZZZZZZZZZ"``, ``"ZZZZZZZZZz"``, ``"ZZZZZZZZzz"``, ``"ZZZZZZZzzz"``, ``"ZZZZZZzzzz"``,`` `` `[`sapply`](https://rdrr.io/r/base/lapply.html)`(`` `` `[`lapply`](https://rdrr.io/r/base/lapply.html)`(``1``:``15``, ``function`` ``(``x``)`` `[`sample`](https://rdrr.io/r/base/sample.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``"Z"``,`[`rep`](https://rdrr.io/r/base/rep.html)`(``"z"``, ``9``)``)``, ``10``)``)``,`` `` ``paste``, collapse``=``""`` `` ``)`` `` ``)``,`` `` XG``=``"CT"`` ``)`` ``#> Writing sample BAM [0.005s]`` ``#> [1] 20`` ``plotMetrics``(``out.bam``, ``as``(``"chrS:1-20"``, ``"GRanges"``)``, title``=``"several overlapping partial epimutations"``)`

![](values_files/figure-html/unnamed-chunk-3-7.png)

` `` ``# simulated long-read sequencing, low methylation`` ``getXM`` ``<-`` ``function`` ``(``p``)`` ``{`[`sample`](https://rdrr.io/r/base/sample.html)`(``x``=`[`c`](https://rdrr.io/r/base/c.html)`(``"z"``, ``"Z"``)``, size``=``1``, prob``=`[`c`](https://rdrr.io/r/base/c.html)`(``p``, ``1``-``p``)``)``}`` ``probs`` ``<-`` ``(`[`sin`](https://rdrr.io/r/base/Trig.html)`(`[`seq`](https://rdrr.io/r/base/seq.html)`(``-``2``*``pi``, ``+``1``*``pi``, by ``=`` ``pi``/``25``)``)``+``2``)``/``3`` `[`simulateBam`](../reference/simulateBam.md)`(`` `` output.bam.file``=``out.bam``,`` `` pos``=``1``:``10``,`` `` XM``=`[`sapply`](https://rdrr.io/r/base/lapply.html)`(``1``:``10``, ``function`` ``(``i``)`` ``{`[`paste`](https://rdrr.io/r/base/paste.html)`(`[`sapply`](https://rdrr.io/r/base/lapply.html)`(``probs``, ``getXM``)``, collapse``=``""``)``}``)``,`` `` XG``=``"CT"`` ``)`` ``#> Writing sample BAM [0.014s]`` ``#> [1] 10`` ``plotMetrics``(``out.bam``, ``as``(``"chrS:1-1000"``, ``"GRanges"``)``, title``=``"simulated long-read sequencing, low methylation"``)`

![](values_files/figure-html/unnamed-chunk-3-8.png)

` ``# simulated long-read sequencing, high methylation`` `[`simulateBam`](../reference/simulateBam.md)`(`` `` output.bam.file``=``out.bam``,`` `` pos``=``1``:``10``,`` `` XM``=`[`sapply`](https://rdrr.io/r/base/lapply.html)`(``1``:``10``, ``function`` ``(``i``)`` ``{`[`paste`](https://rdrr.io/r/base/paste.html)`(`[`sapply`](https://rdrr.io/r/base/lapply.html)`(``1``-``probs``, ``getXM``)``, collapse``=``""``)``}``)``,`` `` XG``=``"CT"`` ``)`` ``#> Writing sample BAM [0.012s]`` ``#> [1] 10`` ``plotMetrics``(``out.bam``, ``as``(``"chrS:1-1000"``, ``"GRanges"``)``, title``=``"simulated long-read sequencing, high methylation"``)`

![](values_files/figure-html/unnamed-chunk-3-9.png)

` ``# amplicon 0%`` ``plotMetrics``(`` `` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(``"extdata"``, ``"amplicon000meth.bam"``, package``=``"epialleleR"``)``,`` `` ``as``(``"chr17:43124861-43126026"``, ``"GRanges"``)``, title``=``"amplicon, 0%"`` ``)`

![](values_files/figure-html/unnamed-chunk-3-10.png)

` ``# amplicon 10%`` ``plotMetrics``(`` `` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(``"extdata"``, ``"amplicon010meth.bam"``, package``=``"epialleleR"``)``,`` `` ``as``(``"chr17:43124861-43126026"``, ``"GRanges"``)``, title``=``"amplicon, 10%"`` ``)`

![](values_files/figure-html/unnamed-chunk-3-11.png)

` ``# sample capture, BMP7`` ``plotMetrics``(`` `` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(``"extdata"``, ``"capture.bam"``, package``=``"epialleleR"``)``,`` `` ``as``(``"chr20:57266125-57268185:+"``, ``"GRanges"``)``, title``=``"sample capture, BMP7, + strand"`` ``)`

![](values_files/figure-html/unnamed-chunk-3-12.png)

` ``# sample capture, BMP7`` ``plotMetrics``(`` `` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(``"extdata"``, ``"capture.bam"``, package``=``"epialleleR"``)``,`` `` ``as``(``"chr20:57266125-57268185:-"``, ``"GRanges"``)``, title``=``"sample capture, BMP7, - strand"`` ``)`

![](values_files/figure-html/unnamed-chunk-3-13.png)

` ``# sample capture, RAD51C`` ``plotMetrics``(`` `` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(``"extdata"``, ``"capture.bam"``, package``=``"epialleleR"``)``,`` `` ``as``(``"chr17:58691673-58693108:+"``, ``"GRanges"``)``, title``=``"sample capture, RAD51C, + strand"`` ``)`

![](values_files/figure-html/unnamed-chunk-3-14.png)

` ``# sample capture, RAD51C`` ``plotMetrics``(`` `` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(``"extdata"``, ``"capture.bam"``, package``=``"epialleleR"``)``,`` `` ``as``(``"chr17:58691673-58693108:-"``, ``"GRanges"``)``, title``=``"sample capture, RAD51C, - strand"`` ``)`

![](values_files/figure-html/unnamed-chunk-3-15.png)

` ``# long-read sequencing, BRCA1`` ``plotMetrics``(`` `` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(``"extdata"``, ``"longread.bam"``, package``=``"epialleleR"``)``,`` `` ``as``(``"chr17:43125000-43127000:+"``, ``"GRanges"``)``, title``=``"long-read sequencing, BRCA1, + strand"``, clip.to.targets``=``TRUE`` ``)`

![](values_files/figure-html/unnamed-chunk-3-16.png)

` ``# long-read sequencing, BRCA1`` ``plotMetrics``(`` `` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(``"extdata"``, ``"longread.bam"``, package``=``"epialleleR"``)``,`` `` ``as``(``"chr17:43125000-43127000:-"``, ``"GRanges"``)``, title``=``"long-read sequencing, BRCA1, - strand"``, clip.to.targets``=``TRUE`` ``)`

![](values_files/figure-html/unnamed-chunk-3-17.png)

### Session Info

[`sessionInfo`](https://rdrr.io/r/utils/sessionInfo.html)`(``)`` ``#> R version 4.6.1 (2026-06-24)`` ``#> Platform: x86_64-pc-linux-gnu`` ``#> Running under: Ubuntu 24.04.4 LTS`` ``#> `` ``#> Matrix products: default`` ``#> BLAS: /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 `` ``#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so; LAPACK version 3.12.0`` ``#> `` ``#> locale:`` ``#> [1] LC_CTYPE=en_US.UTF-8 LC_NUMERIC=C LC_TIME=en_US.UTF-8 `` ``#> [4] LC_COLLATE=en_US.UTF-8 LC_MONETARY=en_US.UTF-8 LC_MESSAGES=en_US.UTF-8 `` ``#> [7] LC_PAPER=en_US.UTF-8 LC_NAME=C LC_ADDRESS=C `` ``#> [10] LC_TELEPHONE=C LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C `` ``#> `` ``#> time zone: UTC`` ``#> tzcode source: system (glibc)`` ``#> `` ``#> attached base packages:`` ``#> [1] stats4 stats graphics grDevices utils datasets methods base `` ``#> `` ``#> other attached packages:`` ``#> [1] epialleleR_1.21.2 ggplot2_4.0.3 GenomicRanges_1.65.4 Seqinfo_1.3.2 `` ``#> [5] IRanges_2.47.5 S4Vectors_0.51.9 BiocGenerics_0.59.12 generics_0.1.4 `` ``#> [9] data.table_1.18.6.1 `` ``#> `` ``#> loaded via a namespace (and not attached):`` ``#> [1] gtable_0.3.6 jsonlite_2.0.0 dplyr_1.2.1 compiler_4.6.1 Rcpp_1.1.2 `` ``#> [6] tidyselect_1.2.1 jquerylib_0.1.4 systemfonts_1.3.2 scales_1.4.0 textshaping_1.0.5 `` ``#> [11] yaml_2.3.12 fastmap_1.2.0 R6_2.6.1 labeling_0.4.3 knitr_1.51 `` ``#> [16] htmlwidgets_1.6.4 tibble_3.3.1 desc_1.4.3 pillar_1.11.1 bslib_0.12.0 `` ``#> [21] RColorBrewer_1.1-3 rlang_1.3.0 cachem_1.1.0 xfun_0.60 fs_2.1.0 `` ``#> [26] sass_0.4.10 S7_0.2.2 otel_0.2.0 cli_3.6.6 withr_3.0.3 `` ``#> [31] pkgdown_2.2.1.9000 magrittr_2.0.5 digest_0.6.39 grid_4.6.1 lifecycle_1.0.5 `` ``#> [36] vctrs_0.7.3 evaluate_1.0.5 glue_1.8.1 farver_2.1.2 ragg_1.5.2 `` ``#> [41] rmarkdown_2.32 pkgconfig_2.0.3 tools_4.6.1 htmltools_0.5.9`

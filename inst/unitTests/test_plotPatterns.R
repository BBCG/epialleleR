test_plotPatterns <- function () {
  amplicon.patterns <- extractPatterns(
    bam=system.file("extdata", "amplicon010meth.bam", package="epialleleR"),
    bed=as("chr17:43124895-43125150", "GRanges"),
    extract.context="CX", highlight.positions="43125043"
  )
  
  selected.patterns <- plotPatterns(
    amplicon.patterns, marginal.transform="log10"
  )
  
  RUnit::checkTrue(
    inherits(selected.patterns, "data.table")
  )
  
  capture.patterns <- extractPatterns(
    bam=system.file("extdata", "capture.bam", package="epialleleR"),
    bed=as("chr17:61864583-61864585", "GRanges"),
    extract.context="CX", highlight.positions=61864584, verbose=FALSE
  )
  
  gtable.patterns <- plotPatterns(
    capture.patterns[strand=="-"], genomic.scale="discrete",  marginal="count",
    plot.context="CxG", tag="pattern", npatterns.per.bin=Inf, plot=FALSE,
    verbose=FALSE
  )
  
  gtable.patterns <- plotPatterns(
    capture.patterns[strand=="+"], genomic.scale="discrete",  marginal="count",
    plot.context="CX", tag="beta", npatterns.per.bin=Inf, plot=FALSE,
    verbose=FALSE
  )
  
  RUnit::checkTrue(
    inherits(gtable.patterns, "gtable")
  )
}

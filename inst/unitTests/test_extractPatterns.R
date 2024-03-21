test_extractPatterns <- function () {
  noclip.patterns <- extractPatterns(
    bam=system.file("extdata", "amplicon010meth.bam", package="epialleleR"),
    bed=system.file("extdata", "amplicon.bed", package="epialleleR"),
    bed.row=2, verbose=TRUE
  )
  
  message(paste(unique(noclip.patterns$pattern), collapse=" "))
  RUnit::checkEquals(
    dim(noclip.patterns),
    c(310, 43)
  )
  
  RUnit::checkEquals(
    length(unique(noclip.patterns$pattern)),
    34
  )
  
  RUnit::checkEquals(
    sum(noclip.patterns$nbase),
    4915
  )
  
  RUnit::checkEquals(
    noclip.patterns[beta>0.5, length(unique(pattern))],
    11
  )
  
  RUnit::checkEquals(
    match(c("43125196", "43125214", "43125957", "43126000"), colnames(noclip.patterns)),
    c(8, 9, 42, 43)
  )
  
  RUnit::checkEquals(
    sum(noclip.patterns=="z", na.rm=TRUE),
    4519
  )
  
  RUnit::checkEquals(
    sum(noclip.patterns=="Z", na.rm=TRUE),
    396
  )
  
  clip.patterns <- extractPatterns(
    bam=system.file("extdata", "amplicon010meth.bam", package="epialleleR"),
    bed=system.file("extdata", "amplicon.bed", package="epialleleR"),
    bed.row=2, clip.patterns=TRUE, verbose=TRUE
  )
  
  RUnit::checkEquals(
    dim(clip.patterns),
    c(154, 26)
  )
  
  RUnit::checkEquals(
    length(unique(clip.patterns$pattern)),
    23
  )
  
  RUnit::checkEquals(
    sum(clip.patterns$nbase),
    2186
  )
  
  RUnit::checkEquals(
    clip.patterns[beta>0.5, length(unique(pattern))],
    8
  )
  
  RUnit::checkEquals(
    sum(clip.patterns=="z", na.rm=TRUE),
    2006
  )
  
  RUnit::checkEquals(
    sum(clip.patterns=="Z", na.rm=TRUE),
    180
  )
  
  exact.patterns <- extractPatterns(
    bam=system.file("extdata", "amplicon010meth.bam", package="epialleleR"),
    bed=as("chr17:43124895-43126001", "GRanges"),
    clip.patterns=TRUE, verbose=TRUE
  )
  
  RUnit::checkEquals(
    length(unique(exact.patterns$pattern)),
    55
  )
  
  RUnit::checkEquals(
    match(c("43124894", "43124895", "43126000", "43126001"), colnames(exact.patterns)),
    c(8, NA, 59, NA)
  )
  
  nooffset.patterns <- extractPatterns(
    bam=system.file("extdata", "amplicon010meth.bam", package="epialleleR"),
    bed=as("chr17:43124895-43126001", "GRanges"),
    clip.patterns=TRUE, strand.offset=0, verbose=TRUE
  )
  
  RUnit::checkEquals(
    length(unique(nooffset.patterns$pattern)),
    55
  )
  
  RUnit::checkEquals(
    match(c("43124894", "43124895", "43126000", "43126001"), colnames(nooffset.patterns)),
    c(NA, 8, NA, 59)
  )
  
  cxg.patterns <- extractPatterns(
    bam=system.file("extdata", "amplicon010meth.bam", package="epialleleR"),
    bed=as("chr17:43124895-43126001", "GRanges"),
    extract.context="CxG", clip.patterns=TRUE, verbose=FALSE
  )
  
  RUnit::checkEquals(
    length(unique(cxg.patterns$pattern)),
    101
  )
  
  RUnit::checkEquals(
    match(c("43124895", "43124896", "43125985", "43126001"), colnames(cxg.patterns)),
    c(8, 9, 126, 127)
  )
  
  cx.patterns <- extractPatterns(
    bam=system.file("extdata", "amplicon010meth.bam", package="epialleleR"),
    bed=as("chr17:43124895-43126001", "GRanges"),
    extract.context="CX", clip.patterns=TRUE, verbose=FALSE
  )
  
  RUnit::checkEquals(
    dim(cx.patterns),
    c(394, 269)
  )
  
  RUnit::checkEquals(
    length(unique(cx.patterns$pattern)),
    135
  )
  
  RUnit::checkEquals(
    match(c("43124895", "43124896", "43125985", "43126001"), colnames(cx.patterns)),
    c(8, 9, 264, 269)
  )
  
  RUnit::checkEquals(
    sum(cx.patterns=="h", na.rm=TRUE),
    18944
  )
  
  RUnit::checkEquals(
    sum(cx.patterns=="H", na.rm=TRUE),
    38
  )
  
  RUnit::checkEquals(
    sum(cx.patterns=="x", na.rm=TRUE),
    8801
  )
  
  RUnit::checkEquals(
    sum(cx.patterns=="X", na.rm=TRUE),
    26
  )
  
  RUnit::checkEquals(
    sum(cx.patterns=="z", na.rm=TRUE),
    5853
  )
  
  RUnit::checkEquals(
    sum(cx.patterns=="Z", na.rm=TRUE),
    565
  )
  
  
  
  capture.patterns <- extractPatterns(
    bam=system.file("extdata", "capture.bam", package="epialleleR"),
    bed=as("chr20:57266125-57268185", "GRanges"),
    verbose=FALSE
  )
  
  RUnit::checkEquals(
    length(unique(capture.patterns$pattern)),
    100
  )
  
  RUnit::checkEquals(
    sum(capture.patterns$nbase),
    1293
  )
  
  RUnit::checkEquals(
    c(nrow(capture.patterns[beta>0.5]), nrow(capture.patterns)),
    c(75, 115)
  )
  
  RUnit::checkEquals(
    capture.patterns[beta>0.5, length(unique(pattern))],
    61
  )
  
  snv.patterns <- extractPatterns(
    bam=system.file("extdata", "capture.bam", package="epialleleR"),
    bed=as("chr17:61864583-61864585", "GRanges"),
    highlight.positions=61864584,
    verbose=FALSE
  )
  
  RUnit::checkEquals(
    match(c("61864475","61864486","61864504","61864584","61864855",
            "61864859","61864871"), colnames(snv.patterns)),
    c(8:14)
  )
  
  RUnit::checkEquals(
    length(unique(snv.patterns$pattern)),
    17
  )
  
  RUnit::checkEquals(
    sum(snv.patterns$nbase),
    55
  )
  
  RUnit::checkEquals(
    c(nrow(snv.patterns[beta>0.5]), nrow(snv.patterns)),
    c(16, 24)
  )
  
  RUnit::checkEquals(
    snv.patterns[beta>0.5, length(unique(pattern))],
    12
  )
  
  RUnit::checkEquals(
    sum(snv.patterns=="C", na.rm=TRUE),
    11
  )
  
  RUnit::checkEquals(
    sum(snv.patterns=="T", na.rm=TRUE),
    8
  )
  
  RUnit::checkEquals(
    sum(snv.patterns=="z", na.rm=TRUE),
    18
  )
  
  RUnit::checkEquals(
    sum(snv.patterns=="Z", na.rm=TRUE),
    37
  )
  
  RUnit::checkEquals(
    snv.patterns[, .N, by=.(strand, `61864584`)][order(strand, `61864584`)]$N,
    c(8, 2, 1, 11, 2)
  )
  
  RUnit::checkEquals(
    snv.patterns[, .N, by=.(strand, `61864584`, pattern)][order(strand, `61864584`)]$N,
    c(3, 1, 2, 1, 1, 1, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 2)
  )
  
  
  RUnit::checkEquals(
    extractPatterns(
      bam=system.file("extdata", "capture.bam", package="epialleleR"),
      bed=as("chr17:61864583-61864585", "GRanges"),
      highlight.positions=c(61864584, 61864584, 61864584),
      verbose=FALSE
    ),
    extractPatterns(
      bam=system.file("extdata", "capture.bam", package="epialleleR"),
      bed=as("chr17:61864583-61864585", "GRanges"),
      highlight.positions=c(61864584, 61864586),
      verbose=FALSE
    )
  )
  
  RUnit::checkEquals(
    extractPatterns(
      bam=system.file("extdata", "capture.bam", package="epialleleR"),
      bed=as("chr17:61864583-61864585", "GRanges"),
      verbose=FALSE
    ),
    extractPatterns(
      bam=system.file("extdata", "capture.bam", package="epialleleR"),
      bed=as("chr17:61864583-61864585", "GRanges"), bed.row=c(1:5),
      highlight.positions=c(1, 2, -61864584),
      verbose=FALSE
    )
  )
  
  RUnit::checkEquals(
    extractPatterns(
      bam=system.file("extdata", "capture.bam", package="epialleleR"),
      bed=as("chr17:61864583-61864585", "GRanges"), bed.row=c(2),
      verbose=FALSE
    ),
    data.table::data.table()
  )
}

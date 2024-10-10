test_generateVcfReport <- function () {
  amplicon.report <- generateVcfReport(
    bam=system.file("extdata", "amplicon010meth.bam", package="epialleleR"),
    bed=system.file("extdata", "amplicon.bed", package="epialleleR"),
    vcf=system.file("extdata", "amplicon.vcf.gz", package="epialleleR"),
    vcf.style="NCBI", verbose=FALSE
  )
  
  capture.vcf <- VariantAnnotation::readVcf(
    system.file("extdata", "capture.vcf.gz", package="epialleleR"))
  
  capture.report <- generateVcfReport(
    bam=system.file("extdata", "capture.bam", package="epialleleR"),
    bed=system.file("extdata", "capture.bed", package="epialleleR"),
    vcf=capture.vcf,
    verbose=FALSE
  )
  
  capture.report.nobed <- generateVcfReport(
    bam=system.file("extdata", "capture.bam", package="epialleleR"),
    bed=NULL,
    vcf=capture.vcf,
    verbose=FALSE
  )
  
  nothreshold.report <- generateVcfReport(
    bam=system.file("extdata", "amplicon010meth.bam", package="epialleleR"),
    bed=system.file("extdata", "amplicon.bed", package="epialleleR"),
    vcf=system.file("extdata", "amplicon.vcf.gz", package="epialleleR"),
    vcf.style="NCBI", threshold.reads=FALSE, verbose=TRUE
  )
  
  RUnit::checkEquals(
    dim(amplicon.report),
    c(56,17)
  )
  
  RUnit::checkTrue(
    identical(capture.report, capture.report.nobed)
  )
  
  RUnit::checkEquals(
    dim(capture.report),
    c(26292,17)
  )
  
  RUnit::checkException(
    generateVcfReport(
      bam=system.file("extdata", "amplicon010meth.bam", package="epialleleR"),
      bed=NULL,
      vcf=system.file("extdata", "amplicon.vcf.gz", package="epialleleR"),
      vcf.style="NCBI", threshold.reads=FALSE, verbose=TRUE
    )
  )
  
  RUnit::checkEquals(
    sum(amplicon.report$`FEp+`, na.rm=TRUE),
    40
  )
  
  RUnit::checkEquals(
    sum(amplicon.report$`FEp-`, na.rm=TRUE),
    40.15024191,
    tolerance=1e-08
  )
  
  RUnit::checkEquals(
    sum(amplicon.report$SumRef, na.rm=TRUE),
    5282
  )
  
  RUnit::checkEquals(
    sum(amplicon.report$SumAlt, na.rm=TRUE),
    14
  )
  
  # some extended consistency checks
  RUnit::checkEquals(
    amplicon.report[, .N, by=.(REF,ALT)][order(REF, ALT)]$N,
    c(3, 4, 1, 7, 2, 13, 11, 5, 4, 2, 3, 1)
  )
  RUnit::checkEquals(
    amplicon.report[, sum(`M+Ref`, na.rm=TRUE), by=.(REF,ALT)][order(REF, ALT)]$V1,
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  )
  RUnit::checkEquals(
    amplicon.report[, sum(`U+Ref`, na.rm=TRUE), by=.(REF,ALT)][order(REF, ALT)]$V1,
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  )
  RUnit::checkEquals(
    amplicon.report[, sum(`M-Ref`, na.rm=TRUE), by=.(REF,ALT)][order(REF, ALT)]$V1,
    c(25, 0, 13, 72, 20, 149, 0, 54, 53, 26, 39, 12)
  )
  RUnit::checkEquals(
    amplicon.report[, sum(`U-Ref`, na.rm=TRUE), by=.(REF,ALT)][order(REF, ALT)]$V1,
    c(201, 0, 142, 722, 196, 1617, 0, 555, 534, 285, 427, 140)
  )
  RUnit::checkEquals(
    amplicon.report[, sum(`M+Alt`, na.rm=TRUE), by=.(REF,ALT)][order(REF, ALT)]$V1,
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  )
  RUnit::checkEquals(
    amplicon.report[, sum(`U+Alt`, na.rm=TRUE), by=.(REF,ALT)][order(REF, ALT)]$V1,
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  )
  RUnit::checkEquals(
    amplicon.report[, sum(`M-Alt`, na.rm=TRUE), by=.(REF,ALT)][order(REF, ALT)]$V1,
    c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0)
  )
  RUnit::checkEquals(
    amplicon.report[, sum(`U-Alt`, na.rm=TRUE), by=.(REF,ALT)][order(REF, ALT)]$V1,
    c(0, 0, 0, 3, 1, 4, 0, 1, 1, 0, 1, 2)
  )
  RUnit::checkEquals(
    amplicon.report[, sum(`SumRef`, na.rm=TRUE), by=.(REF,ALT)][order(REF, ALT)]$V1,
    c(226, 0, 155, 794, 216, 1766, 0, 609, 587, 311, 466, 152)
  )
  RUnit::checkEquals(
    amplicon.report[, sum(`SumAlt`, na.rm=TRUE), by=.(REF,ALT)][order(REF, ALT)]$V1,
    c(0, 0, 0, 4, 1, 4, 0, 1, 1, 0, 1, 2)
  )
  RUnit::checkEquals(
    amplicon.report[, sum(as.numeric(range)), by=.(REF,ALT)][order(REF,ALT)]$V1,
    c(129375557, 172502178, 43125927, 301876978, 86251276, 560631993, 474379441, 215627738, 172503074, 86251776, 129377635, 43125475)
  )
  
  RUnit::checkEquals(
    sum(capture.report$`FEp+`, na.rm=TRUE),
    18217
  )
  
  RUnit::checkEquals(
    sum(capture.report$`FEp-`, na.rm=TRUE),
    18138
  )
  
  # more extended consistency checks
  RUnit::checkEquals(
    capture.report[, .N, by=.(REF,ALT)][order(REF, ALT)]$N,
    c(798, 2526, 664, 1633, 1816, 5598, 5628, 1859, 1740, 643, 2477, 910)
  )
  RUnit::checkEquals(
    capture.report[, sum(`M+Ref`, na.rm=TRUE), by=.(REF,ALT)][order(REF, ALT)]$V1,
    c(63, 340, 84, 158, 205, 0, 736, 172, 197, 83, 0, 118)
  )
  RUnit::checkEquals(
    capture.report[, sum(`U+Ref`, na.rm=TRUE), by=.(REF,ALT)][order(REF, ALT)]$V1,
    c(175, 707, 181, 356, 434, 0, 1215, 412, 399, 169, 0, 258)
  )
  RUnit::checkEquals(
    capture.report[, sum(`M-Ref`, na.rm=TRUE), by=.(REF,ALT)][order(REF, ALT)]$V1,
    c(73, 0, 70, 202, 163, 758, 0, 155, 152, 83, 302, 121)
  )
  RUnit::checkEquals(
    capture.report[, sum(`U-Ref`, na.rm=TRUE), by=.(REF,ALT)][order(REF, ALT)]$V1,
    c(181, 0, 184, 384, 462, 1183, 0, 398, 389, 160, 727, 274)
  )
  RUnit::checkEquals(
    capture.report[, sum(`M+Alt`, `U+Alt`, `M-Alt`, `U-Alt`, `SumAlt`, na.rm=TRUE), by=.(REF,ALT)][order(REF, ALT)]$V1,
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  )
  RUnit::checkEquals(
    capture.report[, sum(`SumRef`, na.rm=TRUE), by=.(REF,ALT)][order(REF, ALT)]$V1,
    c(492, 1047, 519, 1100, 1264, 1941, 1951, 1137, 1137, 495, 1029, 771)
  )
  RUnit::checkEquals(
    capture.report[, sum(as.numeric(range)), by=.(REF,ALT)][order(REF,ALT)]$V1,
    c(56493857996, 183266744117, 47064083930, 115940986680, 125427574153, 382891917031,
      395294093302, 129504774869, 118485346283, 44376646682, 176296716567, 65676743408)
  )
  
  RUnit::checkEquals(
    dim(nothreshold.report),
    c(56,17)
  )
  
  RUnit::checkEquals(
    sum(nothreshold.report$`FEp+`, na.rm=TRUE),
    40
  )
  
  RUnit::checkEquals(
    sum(nothreshold.report$`FEp-`, na.rm=TRUE),
    41
  )
  
  generateVcfReport(
    bam=system.file("extdata", "capture.bam", package="epialleleR"),
    bed=system.file("extdata", "capture.bed", package="epialleleR"),
    vcf=capture.vcf, report.file=tempfile(), verbose=TRUE
  )
  
  quality.report <- generateVcfReport(
    bam=system.file("extdata", "amplicon010meth.bam", package="epialleleR"),
    bed=system.file("extdata", "amplicon.bed", package="epialleleR"),
    vcf=system.file("extdata", "amplicon.vcf.gz", package="epialleleR"),
    vcf.style="NCBI", threshold.reads=FALSE, verbose=TRUE,
    min.mapq=30, min.baseq=20
  )
  
  RUnit::checkEquals(
    sum(quality.report$`FEp+`, na.rm=TRUE),
    40
  )
  
  RUnit::checkEquals(
    sum(quality.report$`FEp-`, na.rm=TRUE),
    41
  )
  
  RUnit::checkEquals(
    sum(quality.report$SumRef, na.rm=TRUE),
    5164
  )
  
  RUnit::checkEquals(
    sum(quality.report$SumAlt, na.rm=TRUE),
    4
  )
  
  libPaths <- .libPaths()
  RUnit::checkException({
    unloadNamespace("VariantAnnotation")
    .libPaths(new="")
    amplicon.report <- generateVcfReport(
      bam=system.file("extdata", "amplicon010meth.bam", package="epialleleR"),
      bed=system.file("extdata", "amplicon.bed", package="epialleleR"),
      vcf=system.file("extdata", "amplicon.vcf.gz", package="epialleleR"),
      vcf.style="NCBI", verbose=FALSE
    )
  })
  .libPaths(new=libPaths)
}

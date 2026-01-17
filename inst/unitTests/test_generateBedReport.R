test_generateBedReport <- function () {
  amplicon.bam    <- system.file("extdata", "amplicon010meth.bam", package="epialleleR")
  amplicon.bed    <- system.file("extdata", "amplicon.bed", package="epialleleR")
  amplicon.report <- generateAmpliconReport(bam=amplicon.bam, bed=amplicon.bed, verbose=FALSE)
  
  nofilter.report <- generateAmpliconReport(bam=amplicon.bam, bed=amplicon.bed, filter.reads=FALSE, verbose=TRUE)
  nothreshold.report <- generateAmpliconReport(bam=amplicon.bam, bed=amplicon.bed, threshold.reads=FALSE, verbose=TRUE)
  
  capture.bam    <- system.file("extdata", "capture.bam", package="epialleleR")
  capture.bed    <- system.file("extdata", "capture.bed", package="epialleleR")
  capture.report <- generateCaptureReport(bam=capture.bam, bed=capture.bed, verbose=FALSE)
  
  RUnit::checkEquals(
    dim(amplicon.report),
    c(5,10)
  )
  
  RUnit::checkEquals(
    dim(capture.report),
    c(565,10)
  )
  
  RUnit::checkEquals(
    sum(amplicon.report$`nreads-`),
    434
  )
  
  RUnit::checkEquals(
    sum(amplicon.report[,.(`nreads+`,`nreads-`)]),
    477
  )
  
  RUnit::checkEquals(
    sum(amplicon.report[,.(nfiltered)]),
    23
  )
  
  RUnit::checkEquals(
    amplicon.report$VEF,
    c(0.08387097, 0.11475409836, 0.05376344086, 0.10714285714, 0.16666667),
  )
  
  RUnit::checkEquals(
    sum(nofilter.report$`nreads-`),
    440
  )
  
  RUnit::checkEquals(
    sum(nofilter.report[,.(`nreads+`,`nreads-`)]),
    500
  )
  
  RUnit::checkEquals(
    nofilter.report$VEF,
    c(0.08333333333, 0.11475409836, 0.05376344086, 0.10714285714, 0.16037736),
  )
  
  RUnit::checkEquals(
    sum(capture.report$`nreads-`, na.rm=TRUE),
    1131
  )
  
  RUnit::checkEquals(
    sum(capture.report[,.(`nreads+`,`nreads-`)], na.rm=TRUE),
    2282
  )
  
  RUnit::checkEquals(
    dim(nothreshold.report),
    c(5,10)
  )
  
  RUnit::checkTrue(
    all(is.na(nothreshold.report$VEF))
  )
  
  RUnit::checkTrue(
    all(is.na(nofilter.report$nfiltered))
  )
  
  generateAmpliconReport(bam=amplicon.bam, bed=amplicon.bed, report.file=tempfile())
  
  
  quality.report <- generateAmpliconReport(bam=amplicon.bam, bed=amplicon.bed,
                                           min.mapq=30, min.baseq=20, verbose=TRUE)
  
  RUnit::checkEquals(
    sum(quality.report$`nreads-`),
    430
  )
  
  RUnit::checkEquals(
    sum(quality.report[,.(`nreads+`,`nreads-`)]),
    470
  )
  
  RUnit::checkEquals(
    quality.report$VEF,
    c(0.08333333333, 0.11475409836, 0.05376344086, 0.10714285714, 0.15789474),
  )
  
  allowall.report <- generateAmpliconReport(bam=amplicon.bam, bed=amplicon.bed, filter.reads=FALSE, threshold.reads=FALSE)
  RUnit::checkTrue(
    all(is.na(c(allowall.report$nfiltered, allowall.report$VEF)))
  )
}

test_generateBedReport <- function () {
  amplicon.bam    <- system.file("extdata", "amplicon010meth.bam", package="epialleleR")
  amplicon.bed    <- system.file("extdata", "amplicon.bed", package="epialleleR")
  amplicon.report <- generateAmpliconReport(bam=amplicon.bam, bed=amplicon.bed, verbose=FALSE)
  
  nothreshold.report <- generateAmpliconReport(bam=amplicon.bam, bed=amplicon.bed, threshold.reads=FALSE, verbose=TRUE)
  
  capture.bam    <- system.file("extdata", "capture.bam", package="epialleleR")
  capture.bed    <- system.file("extdata", "capture.bed", package="epialleleR")
  capture.report <- generateCaptureReport(bam=capture.bam, bed=capture.bed, verbose=FALSE)
  
  RUnit::checkEquals(
    dim(amplicon.report),
    c(5,9)
  )
  
  RUnit::checkEquals(
    dim(capture.report),
    c(565,9)
  )
  
  RUnit::checkEquals(
    sum(amplicon.report$`nreads-`),
    440
  )
  
  RUnit::checkEquals(
    sum(amplicon.report[,.(`nreads+`,`nreads-`)]),
    500
  )
  
  RUnit::checkEquals(
    amplicon.report$VEF,
    c(0.08333333333, 0.11475409836, 0.05376344086, 0.10714285714, 0.13207547170),
  )
  
  RUnit::checkEquals(
    sum(capture.report$`nreads-`, na.rm=TRUE),
    1472
  )
  
  RUnit::checkEquals(
    sum(capture.report[,.(`nreads+`,`nreads-`)], na.rm=TRUE),
    2968
  )
  
  RUnit::checkEquals(
    dim(nothreshold.report),
    c(5,9)
  )
  
  RUnit::checkTrue(
    all(is.na(nothreshold.report$VEF))
  )
  
  generateAmpliconReport(bam=amplicon.bam, bed=amplicon.bed, report.file=tempfile())
  
  
  quality.report <- generateAmpliconReport(bam=amplicon.bam, bed=amplicon.bed,
                                           min.mapq=30, min.baseq=20, verbose=TRUE)
  
  RUnit::checkEquals(
    sum(quality.report$`nreads-`),
    434
  )
  
  RUnit::checkEquals(
    sum(quality.report[,.(`nreads+`,`nreads-`)]),
    485
  )
  
  RUnit::checkTrue(
    sum(amplicon.report[1:4]$VEF) == sum(quality.report[1:4]$VEF)
  )
  
  RUnit::checkTrue(
    amplicon.report[5]$VEF != quality.report[5]$VEF
  )
  
  RUnit::checkEquals(
    quality.report$VEF,
    c(0.08333333333, 0.11475409836, 0.05376344086, 0.10714285714, 0.13186813187),
  )
}

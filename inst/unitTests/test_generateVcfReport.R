test_generateVcfReport <- function () {
  amplicon.report <- generateVcfReport(
    bam=system.file("extdata", "amplicon010meth.bam", package="epialleleR"),
    bed=system.file("extdata", "amplicon.bed", package="epialleleR"),
    vcf=system.file("extdata", "capture.vcf.gz", package="epialleleR"),
    verbose=FALSE
  )
  
  capture.report <- generateVcfReport(
    bam=system.file("extdata", "capture.bam", package="epialleleR"),
    bed=system.file("extdata", "capture.bed", package="epialleleR"),
    vcf=system.file("extdata", "capture.vcf.gz", package="epialleleR"),
    verbose=FALSE
  )
  
  RUnit::checkEquals(
    dim(amplicon.report),
    c(70,17)
  )
  
  RUnit::checkEquals(
    dim(capture.report),
    c(26292,17)
  )
  
  RUnit::checkEquals(
    sum(amplicon.report$`FEp+`, na.rm=TRUE),
    50
  )
  
  RUnit::checkEquals(
    sum(amplicon.report$`FEp-`, na.rm=TRUE),
    45.913162480649,
    tolerance=1e-08
  )
  
  RUnit::checkEquals(
    sum(capture.report$`FEp+`, na.rm=TRUE),
    18217
  )
  
  RUnit::checkEquals(
    sum(capture.report$`FEp-`, na.rm=TRUE),
    18138
  )
}

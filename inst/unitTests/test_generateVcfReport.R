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
  
  RUnit::checkEquals(
    dim(capture.report),
    c(26292,17)
  )
  
  RUnit::checkEquals(
    sum(amplicon.report$`FEp+`, na.rm=TRUE),
    40
  )
  
  RUnit::checkEquals(
    sum(amplicon.report$`FEp-`, na.rm=TRUE),
    39.45658124,
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
}

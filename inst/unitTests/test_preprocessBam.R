test_preprocessBam <- function () {
  capture.bam  <- system.file("extdata", "capture.bam", package="epialleleR")
  capture.data <- preprocessBam(capture.bam, verbose=FALSE)
  RUnit::checkEquals(
    dim(capture.data),
    c(2968,5)
  )
  
  RUnit::checkEquals(
    preprocessBam(capture.data, verbose=TRUE),
    capture.data
  )
  
  amplicon.bam  <- system.file("extdata", "amplicon010meth.bam", package="epialleleR")
  amplicon.data <- preprocessBam(amplicon.bam, skip.duplicates=TRUE, verbose=FALSE)
  RUnit::checkEquals(
    dim(amplicon.data),
    c(500,5)
  )
  
  quality.data <- preprocessBam(capture.bam, verbose=FALSE,
                                min.mapq=30, min.baseq=20)
  RUnit::checkEquals(
    dim(quality.data),
    c(2968,5)
  )
  
  RUnit::checkTrue(
    !identical(capture.data, quality.data)
  )
  
  RUnit::checkTrue(
    identical(capture.data[,1:3], quality.data[,1:3])
  )
}

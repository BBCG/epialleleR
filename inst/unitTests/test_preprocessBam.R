test_preprocessBam <- function () {
  capture.bam  <- system.file("extdata", "capture.bam", package="epialleleR")
  capture.data <- preprocessBam(capture.bam, verbose=FALSE)
  RUnit::checkEquals(
    dim(capture.data),
    c(2968,4)
  )
  
  RUnit::checkEquals(
    preprocessBam(capture.data, verbose=TRUE),
    capture.data
  )
  
  amplicon.bam  <- system.file("extdata", "amplicon010meth.bam", package="epialleleR")
  amplicon.data <- preprocessBam(amplicon.bam, skip.duplicates=TRUE, verbose=FALSE)
  RUnit::checkEquals(
    dim(amplicon.data),
    c(500,4)
  )
  
  quality.data <- preprocessBam(capture.bam, verbose=FALSE,
                                min.mapq=30, min.baseq=20, nthreads=0)
  RUnit::checkEquals(
    dim(quality.data),
    c(2968,4)
  )
  
  RUnit::checkTrue(
    !identical(capture.data, quality.data)
  )
  
  RUnit::checkTrue(
    identical(data.table::data.table(capture.data[,1:3]), data.table::data.table(quality.data[,1:3]))
  )
  
  if (require(Rsamtools, quietly=TRUE)) {
    RUnit::checkException(
      preprocessBam(file.path(system.file("extdata", package="Rsamtools"), "ex1.bam"), verbose=FALSE)
    )
    
    RUnit::checkException(
      preprocessBam(file.path(system.file("extdata", package="Rsamtools"), "tiny.bam"), verbose=FALSE)
    )
  }
  
}

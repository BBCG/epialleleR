test_preprocessBam <- function () {
  
  ### good BAMs
  
  capture.bam  <- system.file("extdata", "capture.bam", package="epialleleR")
  capture.data <- preprocessBam(capture.bam, verbose=FALSE)
  RUnit::checkEquals(
    dim(capture.data),
    c(2968,4)
  )
  
  nil <- preprocessBam(capture.data, verbose=TRUE)
  nil <- preprocessBam(capture.data, verbose=TRUE, min.mapq=10)
  nil <- preprocessBam(capture.data, verbose=TRUE, nthreads=2)
  
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
  
  ### test BAMs
  
  # empty
  RUnit::checkException(
    preprocessBam(system.file("extdata", "test", "empty.bam", package="epialleleR"), verbose=TRUE)
  )
  
  # paired, name-sorted, with XM
  RUnit::checkTrue(
    methods::is(
      preprocessBam(system.file("extdata", "test", "dragen-pe-namesort-xg-xm.bam", package="epialleleR"), verbose=TRUE),
      "data.table"
    )
  )
  
  # paired, name-sorted, no XM
  RUnit::checkException(
    preprocessBam(system.file("extdata", "test", "dragen-pe-namesort-xg.bam", package="epialleleR"), verbose=TRUE)
  )
  
  # paired, unsorted, with XM
  RUnit::checkException(
    preprocessBam(system.file("extdata", "test", "dragen-pe-unsort-xg-xm.bam", package="epialleleR"), verbose=TRUE)
  )
  
  # paired, unsorted, no XM
  RUnit::checkException(
    preprocessBam(system.file("extdata", "test", "dragen-pe-unsort-xg.bam", package="epialleleR"), verbose=TRUE)
  )
  
  # single-ended, unsorted, with XM
  RUnit::checkTrue(
    methods::is(
      preprocessBam(system.file("extdata", "test", "dragen-se-unsort-xg-xm.bam", package="epialleleR"), verbose=TRUE, skip.duplicates=TRUE),
      "data.table"
    )
  )
  
  # single-ended, unsorted, no XG but there's YD
  RUnit::checkException(
    preprocessBam(system.file("extdata", "test", "bwameth-se-unsort-yd.bam", package="epialleleR"), verbose=TRUE)
  )
  
  # single-ended, unsorted, no XG but there's ZS
  RUnit::checkException(
    preprocessBam(system.file("extdata", "test", "bsmap-se-unsort-zs.bam", package="epialleleR"), verbose=TRUE)
  )
  
  # single-ended, unsorted, no XM
  RUnit::checkException(
    preprocessBam(system.file("extdata", "test", "dragen-se-unsort-xg.bam", package="epialleleR"), verbose=TRUE)
  )
  
  # paired-ended instead of single-ended
  RUnit::checkException(
    preprocessBam(system.file("extdata", "test", "dragen-pe-namesort-xg-xm.bam", package="epialleleR"), paired=FALSE, verbose=TRUE)
  )
  
  # single-ended instead of paired-ended
  RUnit::checkException(
    preprocessBam(system.file("extdata", "test", "dragen-se-unsort-xg-xm.bam", package="epialleleR"), paired=TRUE, verbose=TRUE)
  )
 
  # internal coverage
  nil <- epialleleR:::rcpp_read_bam_single(system.file("extdata", "amplicon000meth.bam", package="epialleleR"), 5, 5, TRUE, 0, 0, 1)
  nil <- epialleleR:::rcpp_read_bam_single(system.file("extdata", "amplicon010meth.bam", package="epialleleR"), 5, 5, TRUE, 1, 1, 1)
  nil <- epialleleR:::rcpp_read_bam_single(system.file("extdata", "amplicon100meth.bam", package="epialleleR"), 5, 5, TRUE, 2, 2, 1)
  nil <- epialleleR:::rcpp_read_bam_single(system.file("extdata", "capture.bam", package="epialleleR"), 5, 5, TRUE, 4, 4, 1)

}

test_preprocessBam <- function () {
  
  ### good BAMs
  
  capture.bam  <- system.file("extdata", "capture.bam", package="epialleleR")
  RUnit::checkIdentical(
    epialleleR:::.checkBam(capture.bam, TRUE)[c("paired", "sorted", "tagged")],
    list(paired=TRUE, sorted=TRUE, tagged="XM")
  )
  
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
  RUnit::checkIdentical(
    epialleleR:::.checkBam(amplicon.bam, TRUE)[c("paired", "sorted", "tagged")],
    list(paired=TRUE, sorted=TRUE, tagged="XM")
  )
  
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
  RUnit::checkIdentical(
    epialleleR:::.checkBam(system.file("extdata", "test", "dragen-pe-namesort-xg-xm.bam", package="epialleleR"), TRUE)[c("paired", "sorted", "tagged")],
    list(paired=TRUE, sorted=TRUE, tagged="XM")
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
  RUnit::checkIdentical(
    epialleleR:::.checkBam(system.file("extdata", "test", "dragen-se-unsort-xg-xm.bam", package="epialleleR"), TRUE)[c("paired", "sorted", "tagged")],
    list(paired=FALSE, sorted=FALSE, tagged="XM")
  )
  
  # long-read, single-ended, unsorted
  longread.bam <- system.file("extdata", "longread.bam", package="epialleleR")
  RUnit::checkTrue(
    methods::is(
      preprocessBam(longread.bam, verbose=TRUE),
      "data.table"
    )
  )
  RUnit::checkIdentical(
    epialleleR:::.checkBam(longread.bam, TRUE)[c("paired", "sorted", "tagged")],
    list(paired=FALSE, sorted=FALSE, tagged="MM")
  )
  
  # simulated long-read, single-ended, unsorted
  out.bam <- tempfile(pattern="simulated", fileext=".bam")
  simulateBam(
    output.bam.file=out.bam,
    pos=1,
    cigar=c("1X4899M1H"),
    tlen=4900,
    Mm=c("C+m,0,2,0;G-m,0,0,0;"),
    Ml=list(as.integer(c(102,128,153,138,101,96)))
  )
  RUnit::checkIdentical(
    epialleleR:::.checkBam(out.bam, TRUE)[c("paired", "sorted", "tagged")],
    list(paired=FALSE, sorted=FALSE, tagged="MM")
  )
  
  
  bed <- as("chr17:43124900-43125250", "GRanges")
  # bed file usage
  bam <- preprocessBam(bam=amplicon.bam, targets=system.file("extdata", "amplicon.bed", package="epialleleR"))
  RUnit::checkIdentical(
    as.integer(unname(bam[, summary(start)])),
    as.integer(c(43124861, 43125171, 43125270, 43125299, 43125624, 43125633))
  )
  
  # paired-end, targets
  bam <- preprocessBam(bam=amplicon.bam, targets=bed)
  RUnit::checkIdentical(
    as.integer(unname(bam[, summary(start)])),
    as.integer(c(43124861, 43124861, 43125171, 43125023, 43125171, 43125173))
  )
  # paired-end, targets, clipped
  bam <- preprocessBam(bam=amplicon.bam, targets=bed, clip.to.targets=TRUE)
  RUnit::checkIdentical(
    as.integer(unname(bam[, summary(start)])),
    as.integer(c(43124900, 43124900, 43125171, 43125042, 43125171, 43125173))
  )
  # single-end, targets
  bam <- preprocessBam(bam=amplicon.bam, paired=FALSE, override.check=TRUE, targets=bed)
  RUnit::checkIdentical(
    as.integer(unname(bam[, summary(start)])),
    as.integer(c(43124861, 43124861, 43125024, 43125024, 43125171, 43125173))
  )
  # single-end, targets, clipped
  bam <- preprocessBam(bam=amplicon.bam, paired=FALSE, override.check=TRUE, targets=bed, clip.to.targets=TRUE)
  RUnit::checkIdentical(
    as.integer(unname(bam[, summary(start)])),
    as.integer(c(43124900, 43124900, 43125024, 43125036, 43125171, 43125173))
  )
  # long-read, targets
  bam <- preprocessBam(bam=longread.bam, min.mapq=30, min.baseq=20, min.prob=178, targets=bed)
  RUnit::checkIdentical(
    as.integer(unname(bam[, summary(start)])),
    as.integer(c(43115251, 43119249, 43122146, 43121293, 43123662, 43125166))
  )
  
  # long-read, targets, clipped
  bam <- preprocessBam(bam=longread.bam, min.mapq=30, min.baseq=20, min.prob=178, targets=bed, clip.to.targets=TRUE)
  RUnit::checkIdentical(
    as.integer(unname(bam[, summary(start)])),
    as.integer(c(43124900, 43124900, 43124900, 43124921, 43124900, 43125166))
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
  
  # paired-ended instead of single-ended, but override
  preprocessBam(system.file("extdata", "test", "dragen-pe-namesort-xg-xm.bam", package="epialleleR"), paired=FALSE, override.check=TRUE, verbose=TRUE)
  
  # single-ended instead of paired-ended
  RUnit::checkException(
    preprocessBam(system.file("extdata", "test", "dragen-se-unsort-xg-xm.bam", package="epialleleR"), paired=TRUE, verbose=TRUE)
  )
 
  # internal coverage
  nil <- epialleleR:::rcpp_read_bam_single_all(system.file("extdata", "amplicon000meth.bam", package="epialleleR"), data.table::data.table(), 5, 5, 2820, 0, 0, 1)
  nil <- epialleleR:::rcpp_read_bam_single_all(system.file("extdata", "amplicon010meth.bam", package="epialleleR"), data.table::data.table(), 5, 5, 2820, 1, 1, 1)
  nil <- epialleleR:::rcpp_read_bam_single_all(system.file("extdata", "amplicon100meth.bam", package="epialleleR"), data.table::data.table(), 5, 5, 2820, 2, 2, 1)
  nil <- epialleleR:::rcpp_read_bam_single_all(system.file("extdata", "capture.bam", package="epialleleR"), data.table::data.table(), 5, 5, 2820, 4, 4, 1)
  nil <- epialleleR:::rcpp_read_bam_mm_single_all(system.file("extdata", "amplicon100meth.bam", package="epialleleR"), data.table::data.table(), 5, 5, -1, TRUE, 2820, 4, 4, 1)
  nil <- epialleleR:::rcpp_read_bam_mm_single_all(system.file("extdata", "capture.bam", package="epialleleR"), data.table::data.table(), 5, 5, -1, TRUE, 2820, 4, 4, 1)

}

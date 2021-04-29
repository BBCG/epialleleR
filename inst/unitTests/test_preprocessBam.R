test_preprocessBam <- function () {
  RUnit::checkException(
    epialleleR:::.processBam(list("bam"=list()), FALSE)
  )
  
  capture.bam  <- system.file("extdata", "capture.bam", package="epialleleR")
  capture.data <- preprocessBam(capture.bam, verbose=FALSE)
  RUnit::checkEquals(
    dim(capture.data),
    c(2968,7)
  )
  
  RUnit::checkEquals(
    preprocessBam(capture.data, verbose=TRUE),
    capture.data
  )
  
  amplicon.bam  <- system.file("extdata", "amplicon010meth.bam", package="epialleleR")
  amplicon.data <- preprocessBam(amplicon.bam, skip.duplicates=TRUE, verbose=FALSE)
  RUnit::checkEquals(
    dim(amplicon.data),
    c(500,7)
  )
  
  # as.character(
  #   GenomicAlignments::sequenceLayer(
  #     Biostrings::BStringSet("abcdefghijklmnopqrstuvwxyz"),
  #     "2S4M4=2X1D4I1N6M1H1P", D.letter="?", N.letter="?"
  #   )
  # )
  RUnit::checkEquals(
    rcpp_apply_cigar("2S4M4=2X1D4I1N6M1H1P", "abcdefghijklmnopqrstuvwxyz", "?"),
    "cdefghijkl??qrstuvwxyz"
  )
}

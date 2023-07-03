test_preprocessGenome <- function () {
  genome <- preprocessGenome(system.file("extdata", "test", "reference.fasta.gz", package="epialleleR"))
  
  RUnit::checkTrue(
    is.list(genome)
  )
  
  RUnit::checkTrue(
    !is.null(attr(genome, "rseq_xptr"))
  )
  
  RUnit::checkEquals(
    genome$rname,
    c("ChrA", "ChrB", "ChrC")
  )
  
  RUnit::checkEquals(
    genome$rlen,
    rep(4900, 3)
  )
  
  RUnit::checkException(
    preprocessGenome(system.file("extdata", "test", package="epialleleR"), verbose=TRUE)
  )
  
  RUnit::checkException(
    preprocessGenome(system.file("extdata", "test", "dragen-se.CX_report.gz", package="epialleleR"), verbose=TRUE)
  )
}

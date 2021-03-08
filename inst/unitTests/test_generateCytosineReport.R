test_generateCytosineReport <- function () {
  capture.bam <- system.file("extdata", "capture.bam", package="epialleleR")
  cg.report   <- generateCytosineReport(capture.bam, verbose=FALSE)
  cx.report   <- generateCytosineReport(capture.bam, threshold.reads=FALSE,
                                        report.context="CX", verbose=FALSE)
  
  RUnit::checkEquals(
    dim(cg.report),
    c(15413,7)
  )
  
  RUnit::checkEquals(
    dim(cx.report),
    c(97237,7)
  )
  
  RUnit::checkEquals(
    sum(cg.report$meth),
    4974
  )

  RUnit::checkEquals(
    sum(cx.report$meth),
    6051
  )
}

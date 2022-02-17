test_generateCytosineReport <- function () {
  capture.bam <- system.file("extdata", "capture.bam", package="epialleleR")
  cg.report   <- generateCytosineReport(capture.bam, verbose=TRUE)
  cx.report   <- generateCytosineReport(capture.bam, threshold.reads=FALSE,
                                        report.context="CX", verbose=FALSE)
  
  RUnit::checkEquals(
    nrow(cx.report[duplicated(paste(rname,pos))]),
    0
  )
  
  RUnit::checkEquals(
    as.numeric(table(cx.report$context)[c("CHH", "CHG", "CG")]),
    c(58292, 23486, 15408)
  )
  
  RUnit::checkEquals(
    dim(cg.report),
    c(15408,6)
  )
  
  RUnit::checkEquals(
    dim(cx.report),
    c(97186,6)
  )
  
  RUnit::checkEquals(
    sum(cg.report$meth),
    4974
  )

  RUnit::checkEquals(
    sum(cg.report$unmeth),
    15245
  )
  
  RUnit::checkEquals(
    sum(cx.report$meth),
    6051
  )
  
  RUnit::checkEquals(
    sum(cx.report$unmeth),
    125903
  )
  
  generateCytosineReport(capture.bam, report.file=tempfile())
  
  
  cg.quality  <- generateCytosineReport(capture.bam, verbose=TRUE,
                                        min.mapq=30, min.baseq=20)
  cx.quality  <- generateCytosineReport(capture.bam, threshold.reads=FALSE,
                                        min.mapq=30, min.baseq=20,
                                        report.context="CX", verbose=FALSE)
  
  RUnit::checkEquals(
    dim(cg.quality),
    c(15182,6)
  )
  
  RUnit::checkEquals(
    dim(cx.quality),
    c(95990,6)
  )
  
  RUnit::checkEquals(
    as.numeric(table(cx.quality$context)[c("CHH", "CHG", "CG")]),
    c(57565, 23243, 15182)
  )
  
  RUnit::checkEquals(
    sum(cg.quality$meth),
    4829
  )
  
  RUnit::checkEquals(
    sum(cg.quality$unmeth),
    15040
  )
  
  RUnit::checkEquals(
    sum(cx.quality$meth),
    5866
  )
  
  RUnit::checkEquals(
    sum(cx.quality$unmeth),
    124057
  )
  
}

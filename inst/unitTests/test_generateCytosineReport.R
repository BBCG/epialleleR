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
    c(15197,6)
  )
  
  RUnit::checkEquals(
    dim(cx.quality),
    c(96151,6)
  )
  
  RUnit::checkEquals(
    as.numeric(table(cx.quality$context)[c("CHH", "CHG", "CG")]),
    c(57687, 23267, 15197)
  )
  
  RUnit::checkEquals(
    sum(cg.quality$meth),
    4830
  )
  
  RUnit::checkEquals(
    sum(cg.quality$unmeth),
    15062
  )
  
  RUnit::checkEquals(
    sum(cx.quality$meth),
    5873
  )
  
  RUnit::checkEquals(
    sum(cx.quality$unmeth),
    124333
  )
  
  
  ### single-end
  
  cx.single <- generateCytosineReport(
    system.file("extdata", "test", "dragen-se-unsort-xg-xm.bam", package="epialleleR"),
    threshold.reads=FALSE, report.context="CX", verbose=TRUE
  )
  
  RUnit::checkEquals(
    dim(cx.single),
    c(3236, 6)
  )
  
  RUnit::checkEquals(
    as.numeric(table(cx.single$context)[c("CHH", "CHG", "CG")]),
    c(2165, 802, 269)
  )
  
  RUnit::checkEquals(
    c(sum(cx.single$meth), sum(cx.single$unmeth)),
    c(355, 3599)
  )
  
}

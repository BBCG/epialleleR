test_generateMhlReport <- function () {
  capture.bam <- system.file("extdata", "capture.bam", package="epialleleR")
  
  generateMhlReport(capture.bam, report.file=tempfile())
  generateMhlReport(preprocessBam(capture.bam), report.file=tempfile())
  
  RUnit::checkTrue(
    identical(
      generateMhlReport(capture.bam, max.haplotype.window=1)[, lmhl],
      generateCytosineReport(capture.bam, threshold.reads=FALSE)[, meth/(meth+unmeth)]
    )
  )
  
  mhl.report <- generateMhlReport(capture.bam)
  RUnit::checkEquals(
    c(sum(mhl.report$coverage), sum(mhl.report[strand=="+"]$coverage), sum(mhl.report[strand=="-"]$coverage)),
    c(20219, 10188, 10031)
  )
  RUnit::checkEquals(
    c(sum(mhl.report$length), sum(mhl.report$lmhl)),
    c(229119.960, 2666.456)
  )
  RUnit::checkEquals(
    c(sum(mhl.report[strand=="+"]$length), sum(mhl.report[strand=="+"]$lmhl)),
    c(119605.010, 1281.342)
  )
  RUnit::checkEquals(
    c(sum(mhl.report[strand=="-"]$length), sum(mhl.report[strand=="-"]$lmhl)),
    c(109514.950, 1385.114)
  )
  
  # amplicon 10%
  amplicon.bam <- system.file("extdata", "amplicon010meth.bam", package="epialleleR")
  # without filtering
  mhl.report <- generateMhlReport(amplicon.bam, max.outofcontext.beta=1)
  RUnit::checkEquals(
    c(sum(mhl.report$coverage), sum(mhl.report[strand=="+"]$coverage), sum(mhl.report[strand=="-"]$coverage)),
    c(7081, 342, 6739)
  )
  RUnit::checkEquals(
    c(sum(mhl.report$length), sum(mhl.report$lmhl)),
    c(6060.46765, 45.78637)
  )
  RUnit::checkEquals(
    c(sum(mhl.report[strand=="+"]$length), sum(mhl.report[strand=="+"]$lmhl)),
    c(2380.83333, 34.03206)
  )
  RUnit::checkEquals(
    c(sum(mhl.report[strand=="-"]$length), sum(mhl.report[strand=="-"]$lmhl)),
    c(3679.63432, 11.75431)
  )
  # with default filtering
  mhl.report <- generateMhlReport(amplicon.bam)
  RUnit::checkEquals(
    c(sum(mhl.report$coverage), sum(mhl.report[strand=="+"]$coverage), sum(mhl.report[strand=="-"]$coverage)),
    c(7070, 339, 6731)
  )
  RUnit::checkEquals(
    c(sum(mhl.report$length), sum(mhl.report$lmhl)),
    c(6051.54262, 43.53694)
  )
  RUnit::checkEquals(
    c(sum(mhl.report[strand=="+"]$length), sum(mhl.report[strand=="+"]$lmhl)),
    c(2375.83333, 33.78206)
  )
  RUnit::checkEquals(
    c(sum(mhl.report[strand=="-"]$length), sum(mhl.report[strand=="-"]$lmhl)),
    c(3675.709286, 9.754883)
  )
  
  # amplicon 100%
  amplicon.bam <- system.file("extdata", "amplicon100meth.bam", package="epialleleR")
  RUnit::checkEquals(
    generateMhlReport(amplicon.bam, min.haplotype.length=1, max.haplotype.window=1,
                      min.mapq=30, min.baseq=20, max.outofcontext.beta=1)[, lmhl],
    generateCytosineReport(amplicon.bam, threshold.reads=FALSE,
                           min.mapq=30, min.baseq=20)[, meth/(meth+unmeth)],
    tolerance=0.022992 # because in lMHL we skip unconverted, while in CX we retain them
  )
  
  # simulated
  out.bam <- tempfile(pattern="simulated", fileext=".bam")
  simulateBam(
    output.bam.file=out.bam,
    cigar=c("10000M1H"),
    XM=c(
      paste(sample(c("Z",rep("z", 9)), 10000, replace=TRUE), collapse=""),
      paste(sample(c("Z",rep("z", 9)), 10000, replace=TRUE), collapse="")
    ),
    XG=c("CT")
  )
  cg.beta <- generateCytosineReport(out.bam, threshold.reads=FALSE)
  mhl.report <- generateMhlReport(out.bam, max.haplotype.window=1)
  RUnit::checkEquals(
    c(sum(mhl.report$coverage), sum(mhl.report$length)),
    c(20000, 100000000)
  )
  RUnit::checkIdentical(
    mhl.report[, lmhl],
    cg.beta[, meth/(meth+unmeth)]
  )
}
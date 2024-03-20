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
  RUnit::checkEquals(
    mhl.report[, sum(as.numeric(pos)), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(82104553191, 73818593632, 35293697221, 32465034595, 46183023478, 54803314759, 32606254666, 22820115100, 24827795998,
      10383726558, 34841084399, 22567815161, 22726223278, 33559193423, 24561128398, 23638910883, 34247644829, 16883136930,
      20839322928, 13610664250, 31263366884, 39381904158, 30528905907, 18746270326, 7167501192, 3948042625, 19918640447,
      20184722006, 21953464255, 13692504247, 25372406639, 23683620028, 42935757410, 36160367626, 3547455654, 3234415920,
      7227571922, 12641484839, 11600091024, 11551404414, 3864101423, 2918776285, 3055619996, 4495408567, 19824204867, 19973772765)
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
  RUnit::checkEquals(
    mhl.report[, sum(as.numeric(pos)), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(229479224, 2064272448, 979611677, 3798461436, 494492965, 167726117, 857544290, 1682667646, 113876489, 332507107, 507181268,
      1635050747, 592544083, 135038031, 1720494501, 532289282, 3464993418, 1525917932, 132592227, 169331819, 304723674, 1175165405,
      310180944, 328875251, 1495710567, 2440023361, 2917278582, 57131133, 427731869, 250321582, 121335075, 41504681, 131853312, 135775)
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
  RUnit::checkEquals(
    mhl.report[, sum(as.numeric(pos)), by=.(rname, strand, context)][order(rname, strand, context)]$V1,
    c(229479224, 2064272448, 843199400, 3798461436, 494492965, 167726117, 857544290, 1495501170, 113876489, 332507107, 507181268,
      1635050747, 592544083, 135038031, 1720494501, 532289282, 3464993418, 1525917932, 132592227, 304723674, 1175165405, 310180944,
      328875251, 1495710567, 2440023361, 2917278582, 57131133, 427731869, 250321582, 121335075, 41504681, 131853312, 135775)
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
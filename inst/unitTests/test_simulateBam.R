test_simulateBam <- function () {
  
  out.bam <- tempfile(pattern="simulated", fileext=".bam")
  simulateBam(
    output.bam.file=out.bam,
    pos=1:6,
    XM=c("ZZZzzZZZ", "ZZzzzzZZ"),
    XG=c("CT", "AG"),
    qual="ABCDEFGH",
    rname="chrZ",
    rnext="chrZ"
  )
  cg.beta <- generateCytosineReport(out.bam, threshold.reads=FALSE)
  RUnit::checkEquals(
    dim(cg.beta),
    c(24, 6)
  )
  RUnit::checkEquals(
    c(sum(cg.beta$meth), sum(cg.beta$unmeth)),
    c(30, 18)
  )
  
  simulateBam(
    output.bam.file=out.bam,
    XM=c(
      paste(rep("Z", 10), collapse=""),
      sapply(
        lapply(1:999, function (x) sample(c("Z",rep("z", 9)), 10)),
        paste, collapse=""
      )
    ),
    XG=c("CT")
  )
  cg.vef  <- generateCytosineReport(out.bam, threshold.reads=TRUE)
  RUnit::checkEquals(
    c(sum(cg.vef$meth), sum(cg.vef$unmeth)),
    c(10, 9990)
  )
  
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
  
  simulateBam(
    output.bam.file=out.bam,
    qname="q1",
    flag=c(99, 147),
    cigar=c("10000M1H"),
    XM=c(
      paste(sample(c("Z",rep("z", 9)), 10000, replace=TRUE), collapse=""),
      paste(sample(c("Z",rep("z", 9)), 10000, replace=TRUE), collapse="")
    ),
    XG=c("CT")
  )
  cg.beta <- generateCytosineReport(out.bam, threshold.reads=FALSE)
  
  out.call <- tempfile(pattern="simulated", fileext=".bam")
  simulateBam(
    output.bam.file=out.call,
    pos=1,
    cigar=c("1X2449M1H"),
    rname=c("ChrA", "ChrB", "ChrC"),
    tlen=2450,
    XG=c("CT")
  )
  callMethylation(
    input.bam.file=out.call,
    output.bam.file=out.bam,
    genome=preprocessGenome(system.file("extdata", "test", "reference.fasta.gz", package="epialleleR"))
  )
  cg.beta <- generateCytosineReport(out.bam, threshold.reads=FALSE)
  
}

test_simulateBam <- function () {
  
  out.bam <- tempfile(pattern="simulated", fileext=".bam")
  simulateBam(
    output.bam.file=out.bam,
    pos=1:6,
    XG=c("CT", "AG"),
    NM=1:12
  )
  
  simulateBam(
    qname="a",
    flag=2,
    rname="chrQ",
    pos=c(1,3),
    mapq=45,
    cigar="5M",
    rnext="chrQ",
    pnext=c(3,1),
    tlen=8,
    seq=c("CCCC", "TTTTTTTT"),
    qual=c("FFFF", "DDDDDDDD"),
    verbose=FALSE,
    XM=c("zzZZ", "ZZzzZZzz")
  )
  
  simulateBam(
    pos=1,
    "AB"=1:10,
    Zf=list(c(1.1, -3.3, 1e-4)),
    ZC=list(10:20), Zc=list(-10:0),
    ZS=list(240:260), Zs=list(-260:-240),
    ZI=list(65530:65540), Zi=list(-65540:-65530),
    output.bam.file=out.bam
  )
  
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
    cigar=c("1X4899M1H"),
    rname=c("ChrA", "ChrB", "ChrC"),
    tlen=4900,
    XG=c("CT")
  )
  callMethylation(
    input.bam.file=out.call,
    output.bam.file=out.bam,
    genome=preprocessGenome(system.file("extdata", "test", "reference.fasta.gz", package="epialleleR"))
  )
  cg.beta <- generateCytosineReport(out.bam, threshold.reads=FALSE)
  
}

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
  
  cx.trim <- generateCytosineReport(capture.bam, threshold.reads=FALSE,
                                    trim=3, report.context="CX", verbose=FALSE)
  cx.notrim <- generateCytosineReport(capture.bam, threshold.reads=FALSE,
                                    trim=0, report.context="CX", verbose=FALSE)
  RUnit::checkTrue(
    ! identical(cx.trim, cx.notrim)
  )
  
  RUnit::checkTrue(
    identical(
      data.table::merge.data.table(
        cx.trim[, .(rname, strand, pos, context)],
        cx.trim[, .(rname, strand, pos, context)]
      ),
      data.table::merge.data.table(
        cx.trim[,   .(rname, strand, pos, context)],
        cx.notrim[, .(rname, strand, pos, context)]
      )
    )
  )
  
  
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
  
  cx.trim <- generateCytosineReport(
    system.file("extdata", "test", "dragen-se-unsort-xg-xm.bam", package="epialleleR"),
    threshold.reads=FALSE, trim=1, report.context="CX", verbose=FALSE
  )
  cx.notrim <- generateCytosineReport(
    system.file("extdata", "test", "dragen-se-unsort-xg-xm.bam", package="epialleleR"),
    threshold.reads=FALSE, trim=0, report.context="CX", verbose=FALSE
  )
  
  RUnit::checkTrue(
    ! identical(cx.trim, cx.notrim)
  )
  
  RUnit::checkTrue(
    identical(
      data.table::merge.data.table(
        cx.trim[, .(rname, strand, pos, context)],
        cx.trim[, .(rname, strand, pos, context)]
      ),
      data.table::merge.data.table(
        cx.trim[,   .(rname, strand, pos, context)],
        cx.notrim[, .(rname, strand, pos, context)]
      )
    )
  )
  
  
  
  ### long-read
  output.bam <- tempfile(pattern="output-", fileext=".bam")
  
  # C+m and chebi for other mod
  # modified from htslib/test/base_mods/MM-chebi.sam
  simulateBam(
    flag=c(0),
    seq=c("AGCTCTCCAGAGTCGNACGCCATYCGCGCGCCACCA"),
    pos=1,
    Mm=c("C+m,2,2,1,4,1;C+76792,6,7;N+n,15;"),
    Ml=list(as.integer(c(102,128,153,179,161,187,212,169))),
    output.bam.file=output.bam
  )
  cx.report <- generateCytosineReport(output.bam, threshold.reads=FALSE, report.context="CX")
  RUnit::checkEquals(
    cx.report[, .(strand, pos, context, meth, unmeth)],
    data.table::data.table(
      strand=factor("+", levels=c("+","-")),
      pos=as.integer(c(3, 5, 7, 8, 14, 18, 20, 21, 25, 27, 29, 31, 32, 34, 35)),
      context=factor(c(2, 2, 2, 6, 7, 7, 2, 2, 7, 7, 7, 2, 2, 2, 2), levels=1:7,
                     labels=c("NA1", "CHH", "NA3", "NA4", "NA5", "CHG", "CG")),
      meth=  as.integer(c(0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0)),
      unmeth=as.integer(c(1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1))
    )
  )
  cx.report <- generateCytosineReport(output.bam, threshold.reads=FALSE, report.context="CX",
                                      min.prob=160, highest.prob=FALSE)
  RUnit::checkEquals(
    cx.report[, .(strand, pos, context, meth, unmeth)],
    data.table::data.table(
      strand=factor("+", levels=c("+","-")),
      pos=as.integer(c(3, 5, 7, 8, 14, 18, 20, 21, 25, 27, 29, 31, 32, 34, 35)),
      context=factor(c(2, 2, 2, 6, 7, 7, 2, 2, 7, 7, 7, 2, 2, 2, 2), levels=1:7,
                     labels=c("NA1", "CHH", "NA3", "NA4", "NA5", "CHG", "CG")),
      meth=  as.integer(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1)),
      unmeth=as.integer(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0))
    )
  )
  
  # C+m and chebi for C+m
  # modified from htslib/test/base_mods/MM-chebi.sam
  simulateBam(
    flag=c(0),
    seq=c("AGCTCTCCAGAGTCGNACGCCATYCGCGCGCCACCA"),
    pos=1,
    Mm=c("C+m,2,2,1,4,1;C+27551,6,7;N+n,15;"),
    Ml=list(as.integer(c(102,128,153,179,161,187,212,169))),
    output.bam.file=output.bam
  )
  cx.report <- generateCytosineReport(output.bam, threshold.reads=FALSE, report.context="CX")
  RUnit::checkEquals(
    cx.report[, .(strand, pos, context, meth, unmeth)],
    data.table::data.table(
      strand=factor("+", levels=c("+","-")),
      pos=as.integer(c(3, 5, 7, 8, 14, 18, 20, 21, 25, 27, 29, 31, 32, 34, 35)),
      context=factor(c(2, 2, 2, 6, 7, 7, 2, 2, 7, 7, 7, 2, 2, 2, 2), levels=1:7,
                     labels=c("NA1", "CHH", "NA3", "NA4", "NA5", "CHG", "CG")),
      meth=  as.integer(c(0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1)),
      unmeth=as.integer(c(1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0))
    )
  )
  cx.report <- generateCytosineReport(output.bam, threshold.reads=FALSE, report.context="CX",
                                      min.prob=160)
  RUnit::checkEquals(
    cx.report[, .(strand, pos, context, meth, unmeth)],
    data.table::data.table(
      strand=factor("+", levels=c("+","-")),
      pos=as.integer(c(3, 5, 7, 8, 14, 18, 20, 21, 25, 27, 29, 31, 32, 34, 35)),
      context=factor(c(2, 2, 2, 6, 7, 7, 2, 2, 7, 7, 7, 2, 2, 2, 2), levels=1:7,
                     labels=c("NA1", "CHH", "NA3", "NA4", "NA5", "CHG", "CG")),
      meth=  as.integer(c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1)),
      unmeth=as.integer(c(1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0))
    )
  )
  
  # modification on both strands
  # modified from htslib/test/base_mods/MM-double.sam
  simulateBam(
    flag=c(0),
    seq=c("AGGATCTCTAGCGGATCGGCGGGGGATATGCCATAT"),
    pos=1,
    Mm=c("C+m,1,3,0;G-m,0,2,0,4;G+o,4;"),
    Ml=list(as.integer(c(128,153,179,115,141,166,192,102))),
    output.bam.file=output.bam
  )
  cx.report <- generateCytosineReport(output.bam, threshold.reads=FALSE, report.context="CX")
  RUnit::checkEquals(
    cx.report[meth>0, .(strand, pos, context)],
    data.table::data.table(
      strand=factor(c("-", "+", "-", "-", "-", "+", "+"), levels=c("+","-")),
      pos=as.integer(c(2, 8, 13, 14, 23, 31, 32)),
      context=factor(c(2, 2, 7, 6, 2, 2, 2), levels=1:7,
                     labels=c("NA1", "CHH", "NA3", "NA4", "NA5", "CHG", "CG"))
    )
  )
  
  #modification on + and -, and on the sequenced and opposite strands
  # modified from htslib/test/base_mods/MM-pileup.sam
  simulateBam(
    flag=c(0, 16),
    seq=c("AGCTCTCCAGAGTCGNACGCCATYCGCGCGCCACCA"),
    pos=1,
    Mm=c("C+m,2,2,1,4,1;C+h,6,7;N+n,15,2;",
         "G-m,0,1,4,1,2;G-h,0,7;N-n,17,2;"),
    Ml=list(as.integer(c(128,153,179,204,230,159,6,215,240)),
            as.integer(c(230,204,179,153,128,6,159,240,215))),
    output.bam.file=output.bam
  )
  cx.report <- generateCytosineReport(output.bam, threshold.reads=FALSE, report.context="CX")
  RUnit::checkEquals(
    unname(unlist(cx.report[strand=="-", .(sum(meth), sum(unmeth))])),
    c(0,8)
  )
  RUnit::checkEquals(
    cx.report[strand=="+" & meth>=1, .(pos, context)],
    data.table::data.table(
      pos=as.integer(c(7, 18, 21, 32, 35)),
      context=factor(c(2, 7, 2, 2, 2), levels=1:7,
                     labels=c("NA1", "CHH", "NA3", "NA4", "NA5", "CHG", "CG"))
    )
  )
  
  # modification on '+' and on the sequenced and opposite strands
  # modified from htslib/test/base_mods/MM-orient.sam
  simulateBam(
    flag=c(0),
    seq=c("AGGATCTCTAGCGGATCGGCGGGGGATATGCCATAT"),
    pos=1,
    Mm=c("C+m,2,0,0;G-m,3,1,1;"),
    Ml=list(as.integer(c(128,153,179,128,153,179))),
    output.bam.file=output.bam
  )
  cx.report <- generateCytosineReport(output.bam, threshold.reads=FALSE, report.context="CX")
  RUnit::checkEquals(
    dim(cx.report),
    c(20,6)
  )
  RUnit::checkEquals(
    cx.report[context=="CG", .(strand, pos, meth, unmeth)],
    data.table::data.table(
      strand=factor(c("+", "-", "+", "-", "+", "-"), levels=c("+","-")),
      pos=as.integer(c(12, 13, 17, 18, 20, 21)),
      meth=  as.integer(rep.int(1, 6)),
      unmeth=as.integer(rep.int(0, 6))
    )
  )
  
  # modification on '-' and on the sequenced and opposite strands
  # modified from htslib/test/base_mods/MM-orient.sam
  simulateBam(
    flag=16,
    seq=c("AGGATCTCTAGCGGATCGGCGGGGGATATGCCATAT"),
    pos=1,
    Mm=c("C+m,5,1,1;G-m,2,0,0;"),
    Ml=list(as.integer(c(128,153,179,128,153,179))),
    output.bam.file=output.bam
  )
  cx.report <- generateCytosineReport(output.bam, threshold.reads=FALSE, report.context="CX")
  RUnit::checkEquals(
    dim(cx.report),
    c(20,6)
  )
  RUnit::checkEquals(
    cx.report[context=="CG", .(strand, pos, meth, unmeth)],
    data.table::data.table(
      strand=factor(c("+", "-", "+", "-", "+", "-"), levels=c("+","-")),
      pos=as.integer(c(12, 13, 17, 18, 20, 21)),
      meth=  as.integer(rep.int(1, 6)),
      unmeth=as.integer(rep.int(0, 6))
    )
  )
  
  # simulateBam(
  #   flag=c(0, 16),
  #   seq=c("AGGATCTCTAGCGGATCGGCGGGGGATATGCCATAT",
  #         "ATATGGCATATCCCCCGCCGATCCGCTAGAGATCCT"),
  #   pos=1,
  #   Mm=c("C+m,1,3,0;", "C+m,1,3,0;",
  #        "G-m,0,0,4,3;", "G-m,0,0,4,3;"),
  #   Ml=list(as.integer(c(128,153,179)), as.integer(c(128,153,179)),
  #           as.integer(c(115,141,166,192)), as.integer(c(115,141,166,192))),
  #   output.bam.file=output.bam
  # )
}

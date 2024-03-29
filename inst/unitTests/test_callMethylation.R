test_callMethylation <- function () {
  output.bam <- tempfile(pattern="output-", fileext=".bam")
  genome <- preprocessGenome(system.file("extdata", "test", "reference.fasta.gz", package="epialleleR"))
  
  # just calls
  
  RUnit::checkEquals(
    callMethylation(
      input.bam.file=system.file("extdata", "test", "dragen-pe-namesort-xg-xm.bam", package="epialleleR"),
      output.bam.file=output.bam, genome=genome, nthreads=0, verbose=FALSE
    ),
    list(nrecs=200, ncalled=0)
  )
  RUnit::checkIdentical(
    epialleleR:::.checkBam(bam.file=output.bam, verbose=TRUE)[c("paired", "sorted", "tagged")],
    list(paired=TRUE, sorted=TRUE, tagged="XM")
  )
  
  RUnit::checkEquals(
    callMethylation(
      input.bam.file=system.file("extdata", "test", "dragen-se-unsort-xg.bam", package="epialleleR"),
      output.bam.file=output.bam, genome=genome
    ),
    list(nrecs=100, ncalled=100)
  )
  RUnit::checkIdentical(
    epialleleR:::.checkBam(bam.file=output.bam, verbose=TRUE)[c("paired", "sorted", "tagged")],
    list(paired=FALSE, sorted=FALSE, tagged="XM")
  )
  
  RUnit::checkEquals(
    callMethylation(
      input.bam.file=system.file("extdata", "test", "bwameth-pe-namesort-yd.bam", package="epialleleR"),
      output.bam.file=output.bam, genome=genome, nthreads=0, verbose=FALSE
    ),
    list(nrecs=200, ncalled=170)
  )
  RUnit::checkIdentical(
    epialleleR:::.checkBam(bam.file=output.bam, verbose=TRUE)[c("paired", "sorted", "tagged")],
    list(paired=TRUE, sorted=TRUE, tagged="XM")
  )
  
  RUnit::checkEquals(
    callMethylation(
      input.bam.file=system.file("extdata", "test", "bwameth-se-unsort-yd.bam", package="epialleleR"),
      output.bam.file=output.bam, genome=genome
    ),
    list(nrecs=100, ncalled=73)
  )
  RUnit::checkIdentical(
    epialleleR:::.checkBam(bam.file=output.bam, verbose=TRUE)[c("paired", "sorted", "tagged")],
    list(paired=FALSE, sorted=FALSE, tagged="XM")
  )
  
  RUnit::checkEquals(
    callMethylation(
      input.bam.file=system.file("extdata", "test", "bsmap-pe-namesort-zs.bam", package="epialleleR"),
      output.bam.file=output.bam, genome=genome, nthreads=0, verbose=FALSE
    ),
    list(nrecs=200, ncalled=200)
  )
  RUnit::checkIdentical(
    epialleleR:::.checkBam(bam.file=output.bam, verbose=TRUE)[c("paired", "sorted", "tagged")],
    list(paired=TRUE, sorted=TRUE, tagged="XM")
  )
  
  RUnit::checkEquals(
    callMethylation(
      input.bam.file=system.file("extdata", "test", "bsmap-se-unsort-zs.bam", package="epialleleR"),
      output.bam.file=output.bam, genome=genome
    ),
    list(nrecs=100, ncalled=100)
  )
  RUnit::checkIdentical(
    epialleleR:::.checkBam(bam.file=output.bam, verbose=TRUE)[c("paired", "sorted", "tagged")],
    list(paired=FALSE, sorted=FALSE, tagged="XM")
  )
  
  # calls - errors
  
  RUnit::checkException(
    callMethylation(
      input.bam.file=system.file("extdata", "test", "empty.bam", package="epialleleR"),
      output.bam.file=output.bam, genome=genome, nthreads=0, verbose=FALSE
    )
  )
  
  RUnit::checkException(
    callMethylation(
      input.bam.file=system.file("extdata", "amplicon000meth.bam", package="epialleleR"),
      output.bam.file=output.bam, genome=genome
    )
  )
  
  RUnit::checkException(
    callMethylation(
      input.bam.file=system.file("extdata", "test", "dragen-se-unsort-xg.bam", package="epialleleR"),
      output.bam.file="", genome=genome, nthreads=0, verbose=FALSE
    )
  )
  
  RUnit::checkException(
    callMethylation(
      input.bam.file=system.file("extdata", "test", "bwameth-se-unsort.bam", package="epialleleR"),
      output.bam.file=output.bam, genome=genome, nthreads=0, verbose=FALSE
    )
  )
  
  RUnit::checkException(
    callMethylation(
      input.bam.file=system.file("extdata", "test", "bwameth-se-unsort.bam", package="epialleleR"),
      output.bam.file="", genome=genome, nthreads=0, verbose=FALSE
    )
  )
  
  RUnit::checkException(
    callMethylation(
      input.bam.file="", output.bam.file=output.bam, genome=genome
    )
  )
  
  # calls - compare
  
  callMethylation(
    input.bam.file=system.file("extdata", "test", "dragen-pe-namesort-xg-xm.bam", package="epialleleR"),
    output.bam.file=output.bam, genome=genome, nthreads=1, verbose=FALSE
  )
  cx.ref  <- generateCytosineReport(system.file("extdata", "test", "dragen-pe-namesort-xg-xm.bam", package="epialleleR"),
                                    threshold.reads=FALSE, report.context="CX") 
  cx.call <- generateCytosineReport(output.bam, threshold.reads=FALSE, report.context="CX") 
  RUnit::checkTrue(
    identical(cx.ref, cx.call)
  )
  
  callMethylation(
    input.bam.file=system.file("extdata", "test", "dragen-pe-namesort-xg.bam", package="epialleleR"),
    output.bam.file=output.bam, genome=genome, nthreads=1, verbose=FALSE
  )
  cx.ref  <- generateCytosineReport(system.file("extdata", "test", "dragen-pe-namesort-xg-xm.bam", package="epialleleR"),
                                    threshold.reads=FALSE, report.context="CX") 
  cx.call <- generateCytosineReport(output.bam, threshold.reads=FALSE, report.context="CX") 
  RUnit::checkTrue(
    identical(cx.ref, cx.call)
  )
  
  callMethylation(
    input.bam.file=system.file("extdata", "test", "dragen-se-unsort-xg.bam", package="epialleleR"),
    output.bam.file=output.bam, genome=genome, nthreads=1, verbose=FALSE
  )
  cx.ref  <- generateCytosineReport(system.file("extdata", "test", "dragen-se-unsort-xg-xm.bam", package="epialleleR"),
                                    threshold.reads=FALSE, report.context="CX") 
  cx.call <- generateCytosineReport(output.bam, threshold.reads=FALSE, report.context="CX") 
  RUnit::checkTrue(
    identical(cx.ref, cx.call)
  )
  
  # bwa-meth and DRAGEN: neither SE nor PE are identical
  callMethylation(
    input.bam.file=system.file("extdata", "test", "bwameth-se-unsort-yd.bam", package="epialleleR"),
    output.bam.file=output.bam, genome=genome, nthreads=1, verbose=FALSE
  )
  cx.ref  <- generateCytosineReport(system.file("extdata", "test", "dragen-se-unsort-xg-xm.bam", package="epialleleR"),
                                    threshold.reads=FALSE, report.context="CX") 
  cx.call <- generateCytosineReport(output.bam, threshold.reads=FALSE, report.context="CX") 
  RUnit::checkTrue(
    ! identical(cx.ref, cx.call)
  )
  
  callMethylation(
    input.bam.file=system.file("extdata", "test", "bwameth-pe-namesort-yd.bam", package="epialleleR"),
    output.bam.file=output.bam, genome=genome, nthreads=1, verbose=FALSE
  )
  cx.ref  <- generateCytosineReport(system.file("extdata", "test", "dragen-pe-namesort-xg-xm.bam", package="epialleleR"),
                                    threshold.reads=FALSE, report.context="CX") 
  cx.call <- generateCytosineReport(output.bam, threshold.reads=FALSE, report.context="CX") 
  RUnit::checkTrue(
    ! identical(cx.ref, cx.call)
  )
  
  # BSMAP and DRAGEN: SE is identical, PE is not
  callMethylation(
    input.bam.file=system.file("extdata", "test", "bsmap-se-unsort-zs.bam", package="epialleleR"),
    output.bam.file=output.bam, genome=genome, nthreads=1, verbose=FALSE
  )
  cx.ref  <- generateCytosineReport(system.file("extdata", "test", "dragen-se-unsort-xg-xm.bam", package="epialleleR"),
                                    threshold.reads=FALSE, report.context="CX") 
  cx.call <- generateCytosineReport(output.bam, threshold.reads=FALSE, report.context="CX") 
  RUnit::checkTrue(
    identical(cx.ref, cx.call)
  )

  callMethylation(
    input.bam.file=system.file("extdata", "test", "bsmap-pe-namesort-zs.bam", package="epialleleR"),
    output.bam.file=output.bam, genome=genome, nthreads=1, verbose=FALSE
  )
  cx.ref  <- generateCytosineReport(system.file("extdata", "test", "dragen-pe-namesort-xg-xm.bam", package="epialleleR"),
                                    threshold.reads=FALSE, report.context="CX") 
  cx.call <- generateCytosineReport(output.bam, threshold.reads=FALSE, report.context="CX") 
  RUnit::checkTrue(
    ! identical(cx.ref, cx.call)
  )
}

test_callMethylation <- function () {
  output.bam <- tempfile(pattern="output-", fileext=".bam")
  genome <- preprocessGenome(system.file("extdata", "test", "reference.fasta.gz", package="epialleleR"))
  
  # just calls
  
  RUnit::checkEquals(
    callMethylation(
      input.bam.file=system.file("extdata", "test", "bwameth-pe-namesort-yc.bam", package="epialleleR"),
      output.bam.file=output.bam, genome=genome, nthreads=0, verbose=FALSE
    ),
    list(nrecs=200, ncalled=170)
  )
  
  RUnit::checkEquals(
    callMethylation(
      input.bam.file=system.file("extdata", "test", "bwameth-se-unsort-yc.bam", package="epialleleR"),
      output.bam.file=output.bam, genome=genome
    ),
    list(nrecs=100, ncalled=73)
  )
  
  RUnit::checkEquals(
    callMethylation(
      input.bam.file=system.file("extdata", "test", "dragen-pe-namesort-xg-xm.bam", package="epialleleR"),
      output.bam.file=output.bam, genome=genome, nthreads=0, verbose=FALSE
    ),
    list(nrecs=200, ncalled=0)
  )
  
  RUnit::checkEquals(
    callMethylation(
      input.bam.file=system.file("extdata", "test", "dragen-se-unsort-xg.bam", package="epialleleR"),
      output.bam.file=output.bam, genome=genome
    ),
    list(nrecs=100, ncalled=100)
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
  
}

extractPatterns <- function (bam,
                             bed,
                             bed.row=1,
                             zero.based.bed=FALSE,
                             match.min.overlap=1,
                             extract.context=c("CG", "CHG", "CHH", "CxG", "CX"),
                             min.context.freq=0.01,
                             clip.patterns=FALSE,
                             strand.offset=c("CG"=1, "CHG"=2, "CHH"=0,
                                             "CxG"=0, "CX"=0)[extract.context],
                             min.mapq=0,
                             min.baseq=0,
                             skip.duplicates=FALSE,
                             nthreads=1,
                             verbose=TRUE)
{
  bed.row         <- as.integer(bed.row[1])
  extract.context <- match.arg(extract.context, extract.context)
  strand.offset   <- as.integer(strand.offset[1])
  
  if (!methods::is(bed, "GRanges"))
    bed <- .readBed(bed.file=bed, zero.based.bed=zero.based.bed,
                    verbose=verbose)
  
  bam <- preprocessBam(bam.file=bam, min.mapq=min.mapq, min.baseq=min.baseq,
                       skip.duplicates=skip.duplicates, nthreads=nthreads,
                       verbose=verbose)
  
  patterns <- .getPatterns(
    bam.processed=bam, bed=bed, bed.row=bed.row,
    match.min.overlap=match.min.overlap,
    extract.context=paste0(.context.to.bases[[extract.context]]
                           [c("ctx.meth","ctx.unmeth")], collapse=""),
    min.context.freq=min.context.freq, clip.patterns=clip.patterns,
    strand.offset=strand.offset, verbose=verbose
  )
  
  return(patterns)
}

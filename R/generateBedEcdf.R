#' generateBedEcdf
#'
#' @description
#' This function
#'
#' @details
#' details
#'
#' @param bam BAM file location string OR preprocessed output of
#' \code{\link{preprocessBam}} function.
#' @param bed Browser Extensible Data (BED) file location string OR object of
#' class \code{\linkS4class{GRanges}} holding genomic coordinates for
#' regions of interest. It is used to match sequencing reads to the genomic
#' regions prior to eCDF computation. The seqlevels of BED file/object must
#' match the seqlevels of the BAM file/object used.
#' @param bed.type character string for the type of assay that was used to
#' produce sequencing reads:
#' \itemize{
#'   \item "amplicon" (the default) -- used for amplicon-based next-generation
#'   sequencing when exact coordinates of sequenced fragments are known.
#'   Matching of reads to genomic ranges are then performed by the read's start
#'   or end positions, either of which should be no further than
#'   `match.tolerance` bases away from the start or end position of genomic
#'   ranges given in BED file/\code{\link[GenomicRanges]{GRanges}} object
#'   \item "capture" -- used for capture-based next-generation sequencing when
#'   reads partially overlap with the capture target regions. Read is considered
#'   to match the genomic range when their overlap is more or equal to
#'   `match.min.overlap`. If read matches two or more BED genomic regions, only
#'   the first match is taken (input \code{\link[GenomicRanges]{GRanges}} are
#'   \strong{not} sorted internally)
#' }
#' @param bed.rows integer vector specifying what `bed` regions should be
#' included in the output. If `c(1)` (the default), then function returns eCDF
#' for the first region of `bed`, if NULL - eCDF functions for all `bed`
#' genomic regions as well as for the reads that didn't match any of the regions
#' (last element of the return value; only if there are such reads).
#' @param zero.based.bed param desc
#' @param match.tolerance param desc
#' @param match.min.overlap param desc
#' @param ecdf.context param desc
#' @param min.mapq param desc
#' @param skip.duplicates param desc
#' @param verbose param desc
#' @return return desc
#' @seealso whatelse
#' @examples
#' \dontrun{
#'   sessionInfo()
#' }
#' @export
generateBedEcdf <- function (bam,
                             bed,
                             bed.type=c("amplicon", "capture"),
                             bed.rows=c(1),
                             zero.based.bed=FALSE,
                             match.tolerance=1,
                             match.min.overlap=1,
                             ecdf.context=c("CG", "CHG", "CHH", "CxG", "CX"),
                             min.mapq=0,
                             skip.duplicates=FALSE,
                             verbose=TRUE)
{
  bed.type     <- match.arg(bed.type, bed.type)
  ecdf.context <- match.arg(ecdf.context, ecdf.context)
  
  if (!is(bed, "GRanges"))
    bed <- .readBed(bed.file=bed, zero.based.bed=zero.based.bed,
                    verbose=verbose)

  bam <- preprocessBam(bam.file=bam, min.mapq=min.mapq,
                       skip.duplicates=skip.duplicates, verbose=verbose)
  
  ecdf.list <- .getBedEcdf(
    bam.processed=bam, bed=bed, bed.type=bed.type, bed.rows=bed.rows,
    match.tolerance=match.tolerance, match.min.overlap=match.min.overlap,
    ctx.meth=.context.to.bases[[ecdf.context]][["ctx.meth"]],
    ctx.unmeth=.context.to.bases[[ecdf.context]][["ctx.unmeth"]],
    ooctx.meth=.context.to.bases[[ecdf.context]][["ooctx.meth"]],
    ooctx.unmeth=.context.to.bases[[ecdf.context]][["ooctx.unmeth"]],
    verbose=verbose
  )
  
  return(ecdf.list)
}


#|[c]{}^

# ##############################################################################


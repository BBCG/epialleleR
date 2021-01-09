#' generateBedEcdf
#'
#' @description
#' `generateBedEcdf` desc. Epiallele-aware. GENERIC
#'
#' @details
#' details
#'
#' @param bam param desc
#' @param bed param desc
#' @param bed.type param desc
#' @param bed.rows param desc
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
  
  if (is.character(bam))
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


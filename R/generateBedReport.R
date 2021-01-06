#' generateBedReport
#'
#' @description
#' `generateBedReport` desc. Epiallele-aware. GENERIC
#'
#' @details
#' details
#'
#' @param bam param desc
#' @param bed.file param desc
#' @param report.file param desc
#' @param zero.based.bed param desc
#' @param bed.type param desc
#' @param match.tolerance param desc
#' @param match.min.overlap param desc
#' @param threshold.reads param desc
#' @param threshold.context param desc
#' @param min.context.sites param desc
#' @param min.context.beta param desc
#' @param max.outofcontext.beta param desc
#' @param report.context param desc
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
generateBedReport <- function (bam,
                               bed.file,
                               report.file=NULL,
                               zero.based.bed=FALSE,
                               bed.type=c("amplicon", "capture"),
                               match.tolerance=1,
                               match.min.overlap=1,
                               threshold.reads=TRUE,
                               threshold.context=c("CG", "CHG", "CHH", "CxG", "CX"),
                               min.context.sites=2,
                               min.context.beta=0.5,
                               max.outofcontext.beta=0.1, # double 0 to 1
                               report.context=threshold.context,
                               min.mapq=0,
                               skip.duplicates=FALSE, # http://www.htslib.org/doc/samtools-markdup.html https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/DuplicateMarking_fDG.htm
                               verbose=TRUE)
{
  bed.type          <- match.arg(bed.type, bed.type)
  threshold.context <- match.arg(threshold.context, threshold.context)
  report.context    <- match.arg(report.context, report.context)
  
  bed <- .readBed(bed.file=bed.file, zero.based.bed=zero.based.bed,
                  verbose=verbose)
  
  if (is.character(bam))
    bam <- .readBam(bam.file=bam, min.mapq=min.mapq,
                    skip.duplicates=skip.duplicates, verbose=verbose)
  
  bam.processed <- .processBam(bam=bam, verbose=verbose)
  
  if (threshold.reads) {
    bam.processed$pass <- .thresholdReads(
      bam.processed=bam.processed,
      ctx.meth=.context.to.bases[[threshold.context]][["ctx.meth"]],
      ctx.unmeth=.context.to.bases[[threshold.context]][["ctx.unmeth"]],
      ooctx.meth=.context.to.bases[[threshold.context]][["ooctx.meth"]],
      ooctx.unmeth=.context.to.bases[[threshold.context]][["ooctx.unmeth"]],
      min.context.sites=min.context.sites,
      min.context.beta=min.context.beta,
      max.outofcontext.beta=max.outofcontext.beta,
      verbose=verbose
    )
  } else {
    bam.processed$pass <- TRUE
  }
  
  bed.report <- .getBedReport(
    bam.processed=bam.processed, bed=bed, bed.type=bed.type,
    match.tolerance=match.tolerance, match.min.overlap=match.min.overlap,
    verbose=verbose
  )
  
  if (is.null(report.file))
    return(bed.report)
  else
    .writeReport(report=bed.report, report.file=report.file, gzip=FALSE,
                 verbose=verbose)
}


#|[c]{}^

# ##############################################################################
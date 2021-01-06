#' generateCytosineReport
#'
#' @description
#' `generateCytosineReport` desc. Epiallele-aware
#'
#' @details
#' details
#'
#' @param bam param desc
#' @param report.file param desc
#' @param threshold.reads param desc
#' @param threshold.context param desc
#' @param min.context.sites param desc
#' @param min.context.beta param desc
#' @param max.outofcontext.beta param desc
#' @param report.context param desc
#' @param min.mapq param desc
#' @param skip.duplicates param desc
#' @param gzip param desc
#' @param verbose param desc
#' @return return desc
#' @seealso whatelse
#' @examples
#' \dontrun{
#'   sessionInfo()
#' }
#' @export
generateCytosineReport <- function (bam,
                                    report.file=NULL,
                                    threshold.reads=TRUE,
                                    threshold.context=c("CG", "CHG", "CHH", "CxG", "CX"),
                                    min.context.sites=2,
                                    min.context.beta=0.5,
                                    max.outofcontext.beta=0.1, # double 0 to 1
                                    report.context=threshold.context,
                                    min.mapq=0,
                                    skip.duplicates=FALSE, # http://www.htslib.org/doc/samtools-markdup.html https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/DuplicateMarking_fDG.htm
                                    gzip=FALSE,
                                    verbose=TRUE)
{
  threshold.context <- match.arg(threshold.context, threshold.context)
  report.context    <- match.arg(report.context, report.context)
  
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
  
  cx.report <- .getCytosineReport(
    bam.processed=bam.processed,
    ctx=paste0(.context.to.bases[[report.context]][["ctx.meth"]],
               .context.to.bases[[report.context]][["ctx.unmeth"]]),
    verbose=verbose
  )
  
  if (is.null(report.file))
    return(cx.report)
  else
    .writeReport(report=cx.report, report.file=report.file, gzip=gzip,
                 verbose=verbose)
}


#|[c]{}^

# ##############################################################################

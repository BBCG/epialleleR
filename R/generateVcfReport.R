#' generateVcfReport
#'
#' @description
#' `generateVcfReport` desc. Epiallele-aware. VCF...
#'
#' @details
#' details
#'
#' @param bam param desc
#' @param vcf param desc
#' @param vcf.style param desc
#' @param bed param desc
#' @param zero.based.bed param desc
#' @param threshold.reads param desc
#' @param threshold.context param desc
#' @param min.context.sites param desc
#' @param min.context.beta param desc
#' @param max.outofcontext.beta param desc
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
generateVcfReport <- function (bam,
                               vcf,
                               vcf.style=NULL,
                               bed=NULL,
                               zero.based.bed=FALSE,
                               threshold.reads=TRUE,
                               threshold.context=c("CG", "CHG", "CHH", "CxG", "CX"),
                               min.context.sites=2,
                               min.context.beta=0.5,
                               max.outofcontext.beta=0.1, # double 0 to 1
                               min.mapq=0,
                               skip.duplicates=FALSE,
                               verbose=TRUE)
{
  threshold.context <- match.arg(threshold.context, threshold.context)
  
  if (!is(vcf, "CollapsedVCF")) {
    if (!is.null(bed) & !is(bed, "GRanges"))
      bed <- .readBed(bed.file=bed, zero.based.bed=zero.based.bed, verbose=verbose)
    vcf <- .readVcf(vcf.file=vcf, vcf.style=vcf.style, bed=bed, verbose=verbose)
  } else {
    if (verbose)
      message("Already preprocessed VCF supplied as an input. Options",
              " 'bed' and 'zero.based.bed' will have no effect.")
  }
  
  bam <- preprocessBam(bam.file=bam, min.mapq=min.mapq,
                       skip.duplicates=skip.duplicates, verbose=verbose)
  if (threshold.reads) {
    bam$pass <- .thresholdReads(
      bam.processed=bam,
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
    bam$pass <- TRUE
  }
  
  vcf.report <- .getBaseFreqReport(bam.processed=bam, vcf=vcf, verbose=verbose)
  
  return(vcf.report[,grep("nam|ran|ref|alt|fep",colnames(vcf.report), ignore.case=TRUE), with=FALSE])
}


#|[c]{}^

# ##############################################################################


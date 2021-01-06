#' generateBedReport
#'
#' @description
#' `generateBedReport` desc. Epiallele-aware. GENERIC
#'
#' @details
#' details
#'
#' @param bam param desc
#' @return return desc
#' @seealso whatelse
#' @examples
#'   sessionInfo
# @import parallel
# @importFrom Rsamtools 
#' @export
generateBedReport <- function (bam,
                               bed.file,
                               report.file=NULL,
                               zero.based.bed=TRUE,
                               bed.type=c("amplicon", "capture")
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
  # del afterwards ----------------------------------------------
    # library(Rsamtools)
    # library(GenomicAlignments)
    # library(Biostrings)
    library(dplyr)
    # library(Rcpp)
    report.file <- NULL
    zero.based.bed <- FALSE
    match.tolerance <- 1
    match.min.overlap <- 1
    min.mapq <- 0 #33
    skip.duplicates <- FALSE
    threshold.reads <- TRUE
    threshold.context <- "CG"
    min.context.sites <- 2
    min.context.beta <- 0.5
    max.outofcontext.beta <- 0.01
    report.context <- threshold.context
    verbose <- TRUE
    
    # bed.type <- "amplicon"
    # bam      <- "/scratch/epialleleR/w/WHIP050-D701-D501_S1_L001/WHIP050-D701-D501_S1_L001.bam"
    # bed.file <- "/Home/siv22/oni062/work/scripts/meth-hector/refdata/BRCA1_primary_targets.bed"
    
    bed.type <- "capture"
    bam      <- "/scratch/epialleleR/b/1006-PT2-M_S3_L001/1006-PT2-M_S3_L001_pe.bam"
    bed.file <- "/Home/siv22/oni062/work/scripts/meth-hector/refdata/150109_GRCh38_MET_GENE_v1_EPI_capture_targets.bed"
    
    report.file <- base::sub("(\\.bam)?$", paste0(".",report.context,"_",bed.type,"_report.txt"), bam, ignore.case=TRUE)

    
    Rcpp::sourceCpp("/Home/siv22/oni062/work/packages/epialleleR/src/rcpp_merge_ends.cpp")
    Rcpp::sourceCpp("/Home/siv22/oni062/work/packages/epialleleR/src/rcpp_threshold_reads.cpp")
    # Rcpp::sourceCpp("/Home/siv22/oni062/work/packages/epialleleR/src/rcpp_get_xm_beta.cpp")
    Rcpp::sourceCpp("/Home/siv22/oni062/work/packages/epialleleR/src/rcpp_match_target.cpp")
  # del afterwards ----------------------------------------------
  
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
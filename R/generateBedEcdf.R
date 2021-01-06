#' generateBedEcdf
#'
#' @description
#' `generateBedEcdf` desc. Epiallele-aware. GENERIC
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
generateBedEcdf <- function (bam,
                             bed.file,
                             bed.type=c("amplicon", "capture")
                             bed.rows=c(1),
                             zero.based.bed=TRUE,
                             match.tolerance=1,
                             match.min.overlap=1,
                             ecdf.context=c("CG", "CHG", "CHH", "CxG", "CX"),
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
  bed.rows <- c(NA, 1, 10)
  zero.based.bed <- FALSE
  match.tolerance <- 1
  match.min.overlap <- 1
  min.mapq <- 0 #33
  skip.duplicates <- FALSE
  ecdf.context <- "CG"
  verbose <- TRUE
  
  # bed.type <- "amplicon"
  # bam      <- "/scratch/epialleleR/w/WHIP050-D701-D501_S1_L001/WHIP050-D701-D501_S1_L001.bam"
  # bed.file <- "/Home/siv22/oni062/work/scripts/meth-hector/refdata/BRCA1_primary_targets.bed"
  
  bed.type <- "capture"
  bam      <- "/scratch/epialleleR/b/1006-PT2-M_S3_L001/1006-PT2-M_S3_L001_pe.bam"
  bed.file <- "/Home/siv22/oni062/work/scripts/meth-hector/refdata/150109_GRCh38_MET_GENE_v1_EPI_capture_targets.bed"
  
  Rcpp::sourceCpp("/Home/siv22/oni062/work/packages/epialleleR/src/rcpp_merge_ends.cpp")
  Rcpp::sourceCpp("/Home/siv22/oni062/work/packages/epialleleR/src/rcpp_get_xm_beta.cpp")
  Rcpp::sourceCpp("/Home/siv22/oni062/work/packages/epialleleR/src/rcpp_match_target.cpp")
  # del afterwards ----------------------------------------------
  
  bed.type     <- match.arg(bed.type, bed.type)
  ecdf.context <- match.arg(ecdf.context, ecdf.context)
  
  bed <- .readBed(bed.file=bed.file, zero.based.bed=zero.based.bed,
                  verbose=verbose)
  
  if (is.character(bam))
    bam <- .readBam(bam.file=bam, min.mapq=min.mapq,
                    skip.duplicates=skip.duplicates, verbose=verbose)
  
  bam.processed <- .processBam(bam=bam, verbose=verbose)
  
  ecdf.list <- .getBedEcdf(
    bam.processed=bam.processed, bed=bed, bed.type=bed.type, bed.rows=bed.rows,
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

x <- 1 # out-of-targets reads
x <- 2 # first target in BED
plot(ecdf.list[[x]]$context, col="red", verticals=TRUE, do.points=FALSE, xlim=c(0,1), xlab="per-read beta value", ylab="cumulative density", main=names(ecdf.list[x]))
plot(ecdf.list[[x]]$out.of.context, col="blue", add=TRUE, verticals=TRUE, do.points=FALSE, xlim=c(0,1), xlab="per-read beta value", ylab="cumulative density", main=names(ecdf.list[x]))
abline(v=0.5, col="grey", lty=2)
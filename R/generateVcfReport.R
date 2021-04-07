#' generateVcfReport
#'
#' @description
#' `generateVcfReport` desc. Epiallele-aware. VCF...
#'
#' @details
#' details
#'
#' @param bam BAM file location string OR preprocessed output of
#' \code{\link{preprocessBam}} function
#' @param vcf Variant Call Format (VCF) file location string OR a VCF object
#' returned by \code{\link[VariantAnnotation]{readVcf}}. If VCF object is
#' supplied, its seqlevels must match the seqlevels of the BAM file/object used
#' @param vcf.style string for the seqlevels style of the VCF file, if
#' different from BED file/object. Only has effect when `vcf` param points to
#' the VCF file location and `bed` is not NULL. Possible values:
#' \itemize{
#'   \item NULL (the default) -- seqlevels in BED file/object and VCF file are
#'   the same
#'   \item "NCBI", "UCSC", ... -- valid parameters of
#'   \code{\link[GenomeInfoDb]{seqlevelsStyle}} function
#' }
#' @param bed Browser Extensible Data (BED) file location string OR object of
#' class \code{\link[GenomicRanges]{GRanges}} holding genomic coordinates for
#' regions of interest. It is used to include only the specific genomic ranges
#' when the VCF file is loaded. This option has no effect when VCF object is
#' supplied as a `vcf` param. The seqlevels of BED file/object must match the
#' seqlevels of the BAM file/object used
#' @param report.file file location string to write the VCF report. If NULL
#' (the default) then report is returned as a
#' \code{\link[data.table]{data.table}} object
#' @param zero.based.bed boolean defining if BED coordinates are zero based
#' (default: FALSE)
#' @param threshold.reads boolean defining if sequence reads should be
#' thresholded before counting bases in reference and variant epialleles
#' (default: TRUE). Disabling thresholding is possible but makes no sense in
#' this context, as all the reads will be assigned to the variant epiallele,
#' which will result in Fisher's Exact test p-value of 1 (in columns `FEp+` and
#' `FEP-`)
#' @param threshold.context string defining cytosine methylation context used
#' for thresholdning the reads
#' \itemize{
#'   \item "CG" (the default) -- within-the-context: CpG cytosines (called as
#'   zZ), out-of-context: all the other cytosines (hHxX)
#'   \item "CHG" -- within-the-context: CHG cytosines (xX), out-of-context: hHzZ
#'   \item "CHH" -- within-the-context: CHH cytosines (hH), out-of-context: xXzZ
#'   \item "CxG" -- within-the-context: CG and CHG cytosines (zZxX),
#'   out-of-context: CHH cytosines (hH)
#'   \item "CX" -- all cytosines are considered within-the-context, this
#'   effectively results in no thresholdning
#' }
#' This option has no effect when read thresholding is disabled
#' @param min.context.sites non-negative integer for minimum number of cytosines
#' within the `threshold.context` (default: 2). Reads containing *fewer*
#' within-the-context cytosines are considered completely unmethylated (thus
#' belonging to the reference epiallele). This option has no effect when read
#' thresholding is disabled
#' @param min.context.beta real number in the range [0;1] (default: 0.5). Reads
#' with average beta value for within-the-context cytosines *below* this
#' threshold are considered completely unmethylated (thus belonging to the
#' reference epiallele). This option has no effect when read thresholding is
#' disabled
#' @param max.outofcontext.beta real number in the range [0;1] (default: 0.1).
#' Reads with average beta value for out-of-context cytosines *above* this
#' threshold are considered completely unmethylated (thus belonging to the
#' reference epiallele). This option has no effect when read thresholding is
#' disabled
#' @param min.mapq non-negative integer threshold for minimum read mapping
#' quality (default: 0). Option has no effect if preprocessed BAM data was
#' supplied as an input
#' @param skip.duplicates boolean defining if duplicate aligned reads should be
#' skipped (default: FALSE). Option has no effect if preprocessed BAM data was
#' supplied as an input OR duplicate reads were not marked by alignment software
#' @param gzip boolean to compress the report (default: FALSE)
#' @param verbose boolean to report progress and timings (default: TRUE)
#' @return \code{\link[data.table]{data.table}} object containing VCF report or
#' NULL if report.file was specified. The report columns are:
#' \itemize{
#'   \item name -- variation identifier (e.g. "rs123456789")
#'   \item seqnames -- reference sequence name
#'   \item range -- genomic coordinates of the variation
#'   \item REF -- base at the reference allele
#'   \item ALT -- base at the alternative allele
#'   \item [M|U][+|-][Ref|Alt] -- number of *Ref*erence or *Alt*ernative bases
#'   that were found at this particular position within *M*ethylated or
#'   *U*nmethylated reads that were mapped to "*+*" (forward) or "*-*" (reverse)
#'   DNA strand
#'   \item SumRef -- sum of all *Ref*erence base counts
#'   \item SumAlt -- sum of all *Alt*ernative base counts
#'   \item FEp+ -- Fisher's Exact test p-value for association of a variation
#'   with methylation status of the reads that map to the "*+*" (forward) DNA
#'   strand. Calculated using following contingency table:
#'   
#'   | `M+Ref` | `M+Alt` |
#'   | `U+Ref` | `U+Alt` |
#'   
#'   \item FEp- -- Fisher's Exact test p-value for association of a variation
#'   with methylation status of the reads that map to the "*-*" (reverse) DNA
#'   strand. Calculated using following contingency table:
#'   
#'   | `M-Ref` | `M-Alt` |
#'   | `U-Ref` | `U-Alt` |
#' }
#' @seealso \code{\link{preprocessBam}} for preloading BAM data,
#' \code{\link{generateCytosineReport}} for methylation statistics at the level
#' of individual cytosines, \code{\link{generateBedReport}} for genomic
#' region-based statistics, and `epialleleR` vignettes for the description of
#' usage and sample data
#' @examples
#'   capture.bam <- system.file("extdata", "capture.bam", package="epialleleR")
#'   capture.bed <- system.file("extdata", "capture.bed", package="epialleleR")
#'   capture.vcf <- system.file("extdata", "capture.vcf.gz", package="epialleleR")
#'   vcf.report <- generateVcfReport(bam=capture.bam, bed=capture.bed, vcf=capture.vcf)
#' @export
generateVcfReport <- function (bam,
                               vcf,
                               vcf.style=NULL,
                               bed=NULL,
                               report.file=NULL,
                               zero.based.bed=FALSE,
                               threshold.reads=TRUE,
                               threshold.context=c("CG", "CHG", "CHH", "CxG", "CX"),
                               min.context.sites=2,
                               min.context.beta=0.5,
                               max.outofcontext.beta=0.1, # double 0 to 1
                               min.mapq=0,
                               skip.duplicates=FALSE,
                               gzip=FALSE,
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
  
  vcf.report <- vcf.report[,grep("nam|ran|ref|alt|fep",colnames(vcf.report), ignore.case=TRUE), with=FALSE]
  
  if (is.null(report.file))
    return(vcf.report)
  else
    .writeReport(report=vcf.report, report.file=report.file, gzip=gzip,
                 verbose=verbose)
}


#|[c]{}^

# ##############################################################################


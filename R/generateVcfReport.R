#' generateVcfReport
#'
#' @description
#' This function calculates base frequencies at particular genomic positions and
#' tests their association with the methylation status of the sequencing reads.
#'
#' @details
#' Using BAM reads and sequence variation information as an input,
#' `generateVcfReport` function thresholds the reads according to supplied
#' parameters and calculates the occurence of \strong{Ref}erence and
#' \strong{Alt}ernative bases withing reads, taking into the account DNA strand
#' the read mapped to and average methylation level (epiallele status) of the
#' read.
#' 
#' The information on sequence variation can be supplied as a Variant Call
#' Format (VCF) file location or an object of class VCF, returned by the
#' \code{\link[VariantAnnotation]{readVcf}} function call. As whole-genome VCF
#' files can be extremely large, it is strongly advised to use only relevant
#' subset of their data, prefiltering the VCF object manually before calling
#' `generateVcfReport` or specifying `bed` parameter when `vcf` points to the
#' location of such large VCF file. Please note that all the BAM, BED and VCF
#' files must use the same seqlevels (i.e. chromosome names).
#' 
#' After counting, function checks if certain bases occur more often within
#' reads belonging to certain epialleles using Fisher's Exact test
#' (\code{\link[stats]{fisher.test}}) and reports separate p-values for reads
#' mapped to \strong{"+"} (forward) and \strong{"-"} (reverse) DNA strands.
#' 
#' Please note that currently final report includes only the VCF entries with
#' single-base REF and ALT alleles. This will be fixed in the future.
#'
#' @param bam BAM file location string OR preprocessed output of
#' \code{\link{preprocessBam}} function.
#' @param vcf Variant Call Format (VCF) file location string OR a VCF object
#' returned by \code{\link[VariantAnnotation]{readVcf}} function. If VCF object
#' is supplied, its seqlevels must match the seqlevels of the BAM file/object
#' used.
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
#' class \code{\linkS4class{GRanges}} holding genomic coordinates for
#' regions of interest. It is used to include only the specific genomic ranges
#' when the VCF file is loaded. This option has no effect when VCF object is
#' supplied as a `vcf` param. The seqlevels of BED file/object must match the
#' seqlevels of the BAM file/object used.
#' @param report.file file location string to write the VCF report. If NULL
#' (the default) then report is returned as a
#' \code{\link[data.table]{data.table}} object.
#' @param zero.based.bed boolean defining if BED coordinates are zero based
#' (default: FALSE).
#' @param threshold.reads boolean defining if sequence reads should be
#' thresholded before counting bases in reference and variant epialleles
#' (default: TRUE). Disabling thresholding is possible but makes no sense in
#' this context, as all the reads will be assigned to the variant epiallele,
#' which will result in Fisher's Exact test p-value of 1 (in columns `FEp+` and
#' `FEP-`).
#' @param threshold.context string defining cytosine methylation context used
#' for thresholdning the reads:
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
#' This option has no effect when read thresholding is disabled.
#' @param min.context.sites non-negative integer for minimum number of cytosines
#' within the `threshold.context` (default: 2). Reads containing \strong{fewer}
#' within-the-context cytosines are considered completely unmethylated (thus
#' belonging to the reference epiallele). This option has no effect when read
#' thresholding is disabled.
#' @param min.context.beta real number in the range [0;1] (default: 0.5). Reads
#' with average beta value for within-the-context cytosines \strong{below} this
#' threshold are considered completely unmethylated (thus belonging to the
#' reference epiallele). This option has no effect when read thresholding is
#' disabled.
#' @param max.outofcontext.beta real number in the range [0;1] (default: 0.1).
#' Reads with average beta value for out-of-context cytosines \strong{above}
#' this threshold are considered completely unmethylated (thus belonging to the
#' reference epiallele). This option has no effect when read thresholding is
#' disabled.
#' @param min.mapq non-negative integer threshold for minimum read mapping
#' quality (default: 0). Option has no effect if preprocessed BAM data was
#' supplied as an input.
#' @param skip.duplicates boolean defining if duplicate aligned reads should be
#' skipped (default: FALSE). Option has no effect if preprocessed BAM data was
#' supplied as an input OR duplicate reads were not marked by alignment
#' software.
#' @param gzip boolean to compress the report (default: FALSE).
#' @param verbose boolean to report progress and timings (default: TRUE).
#' @return \code{\link[data.table]{data.table}} object containing VCF report or
#' NULL if report.file was specified. The report columns are:
#' \itemize{
#'   \item name -- variation identifier (e.g. "rs123456789")
#'   \item seqnames -- reference sequence name
#'   \item range -- genomic coordinates of the variation
#'   \item REF -- base at the reference allele
#'   \item ALT -- base at the alternative allele
#'   \item [M|U][+|-][Ref|Alt] -- number of \strong{Ref}erence or
#'   \strong{Alt}ernative bases that were found at this particular position
#'   within \strong{M}ethylated (above threshold) or \strong{U}nmethylated
#'   (below threshold) reads that were mapped to \strong{"+"} (forward) or
#'   \strong{"-"} (reverse) DNA strand
#'   \item SumRef -- sum of all \strong{Ref}erence base counts
#'   \item SumAlt -- sum of all \strong{Alt}ernative base counts
#'   \item FEp+ -- Fisher's Exact test p-value for association of a variation
#'   with methylation status of the reads that map to the \strong{"+"} (forward)
#'   DNA strand. Calculated using following contingency table:
#'   \tabular{ll}{
#'     M+Ref \tab M+Alt \cr
#'     U+Ref \tab U+Alt 
#'   }
#'   
#'   \item FEp- -- Fisher's Exact test p-value for association of a variation
#'   with methylation status of the reads that map to the \strong{"-"} (reverse)
#'   DNA strand. Calculated using following contingency table:
#'   \tabular{ll}{
#'     M-Ref \tab M-Alt \cr
#'     U-Ref \tab U-Alt 
#'   }
#' }
#' @seealso \code{\link{preprocessBam}} for preloading BAM data,
#' \code{\link{generateCytosineReport}} for methylation statistics at the level
#' of individual cytosines, \code{\link{generateBedReport}} for genomic
#' region-based statistics, and `epialleleR` vignettes for the description of
#' usage and sample data.
#' 
#' \code{\link[GenomicRanges]{GRanges}} class for working with genomic ranges,
#' \code{\link[VariantAnnotation]{readVcf}} function for loading VCF data,
#' \code{\link[GenomeInfoDb]{seqlevelsStyle}} function for getting or setting
#' the seqlevels style.
#' @examples
#'   capture.bam <- system.file("extdata", "capture.bam", package="epialleleR")
#'   capture.bed <- system.file("extdata", "capture.bed", package="epialleleR")
#'   capture.vcf <- system.file("extdata", "capture.vcf.gz", package="epialleleR")
#'   
#'   # VCF report
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
                               max.outofcontext.beta=0.1,
                               min.mapq=0,
                               skip.duplicates=FALSE,
                               gzip=FALSE,
                               verbose=TRUE)
{
  threshold.context <- match.arg(threshold.context, threshold.context)
  
  if (!any(is(vcf, "CollapsedVCF"), is(vcf, "ExpandedVCF"))) {
    if (!is.null(bed) & !is(bed, "GRanges"))
      bed <- .readBed(bed.file=bed, zero.based.bed=zero.based.bed, verbose=verbose)
    vcf <- .readVcf(vcf.file=vcf, vcf.style=vcf.style, bed=bed, verbose=verbose)
  } else {
    if (verbose)
      message("Already preprocessed VCF supplied as an input. Options",
              " 'bed' and 'zero.based.bed' will have no effect.")
  }
  if (is(vcf, "CollapsedVCF"))
    vcf <- VariantAnnotation::expand(vcf, row.names=TRUE)
  
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


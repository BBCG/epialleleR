#' @name generateBedReport
#' @aliases generateAmpliconReport
#' @aliases generateCaptureReport
#'
#' @title generateBedReport
#'
#' @description
#' `generateBedReport`, `generateAmpliconReport`, `generateCaptureReport` --
#' these functions match BAM reads to the set of genomic locations and return
#' the fraction of reads with an average methylation level passing an arbitrary
#' threshold.
#'
#' @details
#' Functions report hypermethylated variant epiallele frequencies (VEF) per
#' genomic region of interest using BAM and BED files as input. Reads (for
#' paired-end sequencing alignment files - read pairs as a single
#' entity) are matched to genomic locations by exact
#' coordinates (`generateAmpliconReport` or `generateBedReport` with an option
#' bed.type="amplicon") or minimum overlap (`generateCaptureReport` or
#' `generateBedReport` with an option bed.type="capture") -- the former to be
#' used for amplicon-based NGS data, while the latter -- for the capture-based
#' NGS data. The function's logic is explained below.
#' 
#' Let's suppose we have a BAM file with four reads, all mapped to the "+"
#' strand of chromosome 1, positions 1-16. The genomic range is supplied as a
#' parameter `bed = as("chr1:1-100", "GRanges")`. Assuming the default values
#' for the thresholding parameters (threshold.reads = TRUE,
#' threshold.context = "CG", min.context.sites = 2, min.context.beta = 0.5,
#' max.outofcontext.beta = 0.1), the input and results will look as following:
#' 
#' \tabular{lll}{
#'   methylation string \tab threshold \tab explained \cr
#'   ...Z..x+.h..x..h. \tab below \tab min.context.sites < 2 (only one zZ base) \cr
#'   ...Z..z.h..x..h.  \tab above \tab pass all criteria \cr
#'   ...Z..z.h..X..h.  \tab below \tab max.outofcontext.beta > 0.1 (1XH / 3xXhH = 0.33) \cr
#'   ...Z..z.h..z-.h.  \tab below \tab min.context.beta < 0.5 (1Z / 3zZ = 0.33)
#' }
#' 
#' Only the second read will satisfy all of the thresholding criteria, leading
#' to the following BED report (given that all reads map to chr1:+:1-16):
#' 
#' \tabular{llllllll}{
#'   seqnames \tab start \tab end \tab width \tab strand \tab nreads+ \tab nreads- \tab VEF \cr
#'   chr1 \tab 1 \tab 100 \tab 100 \tab * \tab 4 \tab 0 \tab 0.25
#' }
#' 
#' Please note, that read thresholding by an average methylation level
#' (as explained above) makes little sense for long-read sequencing alignments,
#' as such reads can cover multiple regions with very different DNA methylation
#' properties. Instead, please use \code{\link{extractPatterns}}, limiting
#' pattern output to the region of interest only.
#' 
#' @param bam BAM file location string OR preprocessed output of
#' \code{\link[epialleleR]{preprocessBam}} function. Read more about BAM file
#' requirements and BAM preprocessing at \code{\link{preprocessBam}}.
#' @param bed Browser Extensible Data (BED) file location string OR object of
#' class \code{\link[GenomicRanges]{GRanges}} holding genomic coordinates for
#' regions of interest. The style of seqlevels of BED file/object must be the
#' same as the style of seqlevels of BAM file/object used. The 
#' BED/\code{\link[GenomicRanges]{GRanges}} rows are \strong{not} sorted
#' internally. As of now, the strand information is ignored and reads (read
#' pairs) matching both strands are separately counted and reported.
#' @param report.file file location string to write the BED report. If NULL
#' (the default) then report is returned as a
#' \code{\link[data.table]{data.table}} object.
#' @param zero.based.bed boolean defining if BED coordinates are zero based
#' (default: FALSE).
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
#' @param match.tolerance integer for the largest difference between read's and
#' BED \code{\link[GenomicRanges]{GRanges}} start or end positions during
#' matching of amplicon-based NGS reads (default: 1).
#' @param match.min.overlap integer for the smallest overlap between read's and
#' BED \code{\link[GenomicRanges]{GRanges}} start or end positions during
#' matching of capture-based NGS reads (default: 1). If read matches two or more
#' BED genomic regions, only the first match is taken (input
#' \code{\link[GenomicRanges]{GRanges}} are \strong{not} sorted internally).
#' @param threshold.reads boolean defining if sequence reads should be
#' thresholded before counting reads belonging to variant epialleles (default:
#' TRUE). Disabling thresholding is possible but makes no sense in the context
#' of this function, because
#' all the reads will be assigned to the variant epiallele, which will result
#' in VEF==1 (in such case `NA` VEF values are returned in order to avoid
#' confusion). As thresholding is \strong{not} recommended for long-read
#' sequencing data, this function is \strong{not} recommended for such data
#' either.
#' @param threshold.context string defining cytosine methylation context used
#' for thresholding the reads:
#' \itemize{
#'   \item "CG" (the default) -- within-the-context: CpG cytosines (called as
#'   zZ), out-of-context: all the other cytosines (hHxX)
#'   \item "CHG" -- within-the-context: CHG cytosines (xX), out-of-context: hHzZ
#'   \item "CHH" -- within-the-context: CHH cytosines (hH), out-of-context: xXzZ
#'   \item "CxG" -- within-the-context: CG and CHG cytosines (zZxX),
#'   out-of-context: CHH cytosines (hH)
#'   \item "CX" -- all cytosines are considered within-the-context, this
#'   effectively results in no thresholding
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
#' @param ... other parameters to pass to the
#' \code{\link[epialleleR]{preprocessBam}} function.
#' Options have no effect if preprocessed BAM data was supplied as an input.
#' @param gzip boolean to compress the report (default: FALSE).
#' @param verbose boolean to report progress and timings (default: TRUE).
#' @return \code{\link[data.table]{data.table}} object containing VEF report for
#' BED \code{\link[GenomicRanges]{GRanges}} or NULL if report.file was
#' specified. If BAM file contains reads that would not match to any of BED
#' \code{\link[GenomicRanges]{GRanges}}, the last row in the report will
#' contain information on such reads (with seqnames, start and end equal to NA).
#' The report columns are:
#' \itemize{
#'   \item seqnames -- reference sequence name
#'   \item start -- start of genomic region
#'   \item end -- end of genomic region
#'   \item width -- width of genomic region
#'   \item strand -- strand
#'   \item ... -- other columns that were present in BED or metadata columns of
#'   \code{\link[GenomicRanges]{GRanges}} object 
#'   \item nreads+ -- number of reads (pairs) mapped to the forward ("+") strand
#'   \item nreads- -- number of reads (pairs) mapped to the reverse ("-") strand
#'   \item VEF -- frequency of reads passing the threshold
#' }
#' @seealso \code{\link{preprocessBam}} for preloading BAM data,
#' \code{\link{generateCytosineReport}} for methylation statistics at the level
#' of individual cytosines, \code{\link{generateVcfReport}} for evaluating
#' epiallele-SNV associations, \code{\link{extractPatterns}} for exploring
#' methylation patterns, \code{\link{generateBedEcdf}} for analysing the
#' distribution of per-read beta values, and `epialleleR` vignettes for the
#' description of usage and sample data.
#' 
#' \code{\link[GenomicRanges]{GRanges}} class for working with genomic ranges,
#' \code{\link[GenomeInfoDb]{seqlevelsStyle}} function for getting or setting
#' the seqlevels style.
#' @examples
#'   # amplicon data
#'   amplicon.bam    <- system.file("extdata", "amplicon010meth.bam",
#'                                  package="epialleleR")
#'   amplicon.bed    <- system.file("extdata", "amplicon.bed",
#'                                  package="epialleleR")
#'   amplicon.report <- generateAmpliconReport(bam=amplicon.bam,
#'                                             bed=amplicon.bed)
#'   
#'   # capture NGS
#'   capture.bam    <- system.file("extdata", "capture.bam",
#'                                 package="epialleleR")
#'   capture.bed    <- system.file("extdata", "capture.bed",
#'                                 package="epialleleR")
#'   capture.report <- generateCaptureReport(bam=capture.bam, bed=capture.bed)
#'   
#'   # generateAmpliconReport and generateCaptureReport are just aliases
#'   # of the generateBedReport
#'   bed.report <- generateBedReport(bam=capture.bam, bed=capture.bed,
#'                                   bed.type="capture")
#'   identical(capture.report, bed.report)
#' @rdname generateBedReport
#' @export
generateAmpliconReport <- function (
  bam, bed, report.file=NULL, zero.based.bed=FALSE, match.tolerance=1,
  threshold.reads=TRUE, threshold.context=c("CG", "CHG", "CHH", "CxG", "CX"),
  min.context.sites=2, min.context.beta=0.5, max.outofcontext.beta=0.1,
  ..., gzip=FALSE, verbose=TRUE)
{
  generateBedReport(
    bam=bam, bed=bed, report.file=report.file, zero.based.bed=zero.based.bed,
    bed.type="amplicon", match.tolerance=match.tolerance,
    threshold.reads=threshold.reads, threshold.context=threshold.context,
    min.context.sites=min.context.sites, min.context.beta=min.context.beta,
    max.outofcontext.beta=max.outofcontext.beta, ..., gzip=gzip, verbose=verbose
  )
}
#' @rdname generateBedReport
#' @export
generateCaptureReport <- function (
  bam, bed, report.file=NULL, zero.based.bed=FALSE, match.min.overlap=1,
  threshold.reads=TRUE, threshold.context=c("CG", "CHG", "CHH", "CxG", "CX"),
  min.context.sites=2, min.context.beta=0.5, max.outofcontext.beta=0.1,
  ..., gzip=FALSE, verbose=TRUE)
{
  generateBedReport(
    bam=bam, bed=bed, report.file=report.file, zero.based.bed=zero.based.bed,
    bed.type="capture", match.min.overlap=match.min.overlap,
    threshold.reads=threshold.reads, threshold.context=threshold.context,
    min.context.sites=min.context.sites, min.context.beta=min.context.beta,
    max.outofcontext.beta=max.outofcontext.beta, ..., gzip=gzip, verbose=verbose
  )
}
#' @rdname generateBedReport
#' @export
generateBedReport <- function (bam,
                               bed,
                               report.file=NULL,
                               zero.based.bed=FALSE,
                               bed.type=c("amplicon", "capture"),
                               match.tolerance=1,
                               match.min.overlap=1,
                               threshold.reads=TRUE,
                               threshold.context=c("CG", "CHG", "CHH", "CxG", "CX"),
                               min.context.sites=2,
                               min.context.beta=0.5,
                               max.outofcontext.beta=0.1,
                               ...,
                               gzip=FALSE,
                               verbose=TRUE)
{
  bed.type          <- match.arg(bed.type, bed.type)
  threshold.context <- match.arg(threshold.context, threshold.context)
  
  if (!methods::is(bed, "GRanges"))
    bed <- .readBed(bed.file=bed, zero.based.bed=zero.based.bed,
                    verbose=verbose)
  
  bam <- preprocessBam(bam.file=bam, ..., verbose=verbose)
  
  if (threshold.reads) {
    pass <- .thresholdReads(
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
    pass <- rep(TRUE, nrow(bam))
  }
  
  bed.report <- .getBedReport(
    bam.processed=bam, pass=pass, bed=bed, bed.type=bed.type,
    match.tolerance=match.tolerance, match.min.overlap=match.min.overlap,
    verbose=verbose
  )
  
  if (!threshold.reads) bed.report$VEF <- NA
  
  if (is.null(report.file))
    return(bed.report)
  else
    .writeReport(report=bed.report, report.file=report.file, gzip=gzip,
                 verbose=verbose)
}

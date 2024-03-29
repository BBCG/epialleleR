#' generateCytosineReport
#'
#' @description
#' This function counts methylated and unmethylated DNA bases taking into the
#' account average methylation level of the entire sequence read.
#'
#' @details
#' The function reports cytosine methylation information using BAM file or data
#' as an input. In contrast to the other currently available software, reads
#' (for paired-end sequencing alignment files - read
#' pairs as a single entity) can be thresholded by their average
#' methylation level before counting methylated bases, effectively resulting in
#' hypermethylated variant epiallele frequency (VEF) being reported instead of
#' beta value. The function's logic is explained below.
#' 
#' Let's suppose we have a BAM file with four reads, all mapped to the "+"
#' strand of chromosome 1, positions 1-16. Assuming the default values
#' for the thresholding parameters (threshold.reads = TRUE,
#' threshold.context = "CG", min.context.sites = 2, min.context.beta = 0.5,
#' max.outofcontext.beta = 0.1), the input and results will look as following:
#' 
#' \tabular{llll}{
#'   methylation string \tab threshold \tab explained \tab methylation reported \cr
#'   ...Z..x+.h..x..h. \tab below \tab min.context.sites < 2 (only one zZ base) \tab all cytosines unmethylated \cr
#'   ...Z..z.h..x..h.  \tab above \tab pass all criteria \tab only C4 (Z at position 4) is methylated \cr
#'   ...Z..z.h..X..h.  \tab below \tab max.outofcontext.beta > 0.1 (1XH / 3xXhH = 0.33) \tab all cytosines unmethylated \cr
#'   ...Z..z.h..z-.h.  \tab below \tab min.context.beta < 0.5 (1Z / 3zZ = 0.33) \tab all cytosines unmethylated
#' }
#' 
#' Only the second read will satisfy all of the thresholding criteria, leading
#' to the following CX report (given that all reads map to chr1:+:1-16):
#' 
#' \tabular{llllll}{
#'   rname \tab strand \tab pos \tab context \tab meth \tab unmeth \cr
#'   chr1 \tab + \tab 4 \tab CG \tab 1 \tab 3 \cr
#'   chr1 \tab + \tab 7 \tab CG \tab 0 \tab 3 \cr
#'   chr1 \tab + \tab 9 \tab CHH \tab 0 \tab 4 \cr
#'   chr1 \tab + \tab 12 \tab CHG \tab 0 \tab 3 \cr
#'   chr1 \tab + \tab 15 \tab CHH \tab 0 \tab 4 
#' }
#' 
#' With the thresholding disabled (threshold.reads = FALSE) all methylated bases
#' will retain their status, so the CX report will be similar to the reports
#' produced by other methylation callers (such as Bismark or Illumina DRAGEN Bio
#' IT Platform):
#' 
#' \tabular{llllll}{
#'   rname \tab strand \tab pos \tab context \tab meth \tab unmeth \cr
#'   chr1 \tab + \tab 4 \tab CG \tab 4 \tab 0 \cr
#'   chr1 \tab + \tab 7 \tab CG \tab 0 \tab 3 \cr
#'   chr1 \tab + \tab 9 \tab CHH \tab 0 \tab 4 \cr
#'   chr1 \tab + \tab 12 \tab CHG \tab 1 \tab 2 \cr
#'   chr1 \tab + \tab 15 \tab CHH \tab 0 \tab 4 
#' }
#' 
#' Other notes:
#' 
#' To produce conventional cytosine reports without thresholding by
#' within-context methylation level though
#' minimally affected by incomplete cytosine conversion, run this method with
#' the following parameters: `threshold.reads=TRUE`, `threshold.context="CG"`,
#' `min.context.sites=0`, `min.context.beta=0`, `max.outofcontext.beta=0.1`.
#' All cytosines within reads (read pairs) having more than 10% out-of-context
#' cytosines methylated, will be effectively treated as unmethylated ones.
#' 
#' Methylation string bases in unknown context ("uU") are simply ignored, which,
#' to the best of our knowledge, is consistent with the behaviour of other
#' tools.
#' 
#' In order to mitigate the effect of sequencing errors (leading to rare
#' variations in the methylation context, as in reads 1 and 4 above), the
#' context present in more than 50\% of the reads is assumed to be correct,
#' while all bases at the same position but having other methylation context
#' are simply ignored. This allows reports to be prepared without using the
#' reference genome sequence.
#' 
#' The downside of not using the reference genome sequence is the inability to
#' determine the actual sequence of triplet for every base in the cytosine 
#' report. Therefore this sequence is not reported, and this won't change
#' until such information will be considered as worth adding.
#'
#' Please also note, that read thresholding by an average methylation level
#' (as explained above) makes little sense for long-read sequencing alignments,
#' as such reads can cover multiple regions with very different DNA methylation
#' properties.
#' 
#' @param bam BAM file location string OR preprocessed output of
#' \code{\link[epialleleR]{preprocessBam}} function. Read more about BAM file
#' requirements and BAM preprocessing at \code{\link{preprocessBam}}.
#' @param report.file file location string to write the cytosine report. If NULL
#' (the default) then report is returned as a
#' \code{\link[data.table]{data.table}} object.
#' @param threshold.reads boolean defining if sequence reads (read pairs) should
#' be thresholded before counting methylated cytosines (default: TRUE).
#' Disabling thresholding makes the report virtually indistinguishable from the
#' ones generated by other software, such as Bismark or Illumina DRAGEN Bio IT
#' Platform. Thresholding is \strong{not} recommended for long-read sequencing
#' data.
#' @param threshold.context string defining cytosine methylation context used
#' for thresholding the reads:
#' \itemize{
#'   \item "CG" (the default) --- within-the-context: CpG cytosines (called as
#'   zZ), out-of-context: all the other cytosines (hHxX)
#'   \item "CHG" --- within-the-context: CHG cytosines (xX), out-of-context: hHzZ
#'   \item "CHH" --- within-the-context: CHH cytosines (hH), out-of-context: xXzZ
#'   \item "CxG" --- within-the-context: CG and CHG cytosines (zZxX),
#'   out-of-context: CHH cytosines (hH)
#'   \item "CX" --- all cytosines are considered within-the-context, this
#'   effectively results in no thresholding
#' }
#' This option has no effect when read thresholding is disabled.
#' @param min.context.sites non-negative integer for minimum number of cytosines
#' within the `threshold.context` (default: 2). Reads containing \strong{fewer}
#' within-the-context cytosines are considered completely unmethylated (all C
#' are counted as T). This option has no effect when read thresholding is
#' disabled.
#' @param min.context.beta real number in the range [0;1] (default: 0.5). Reads
#' with average beta value for within-the-context cytosines \strong{below} this
#' threshold are considered completely unmethylated (all C are counted as T).
#' This option has no effect when read thresholding is disabled.
#' @param max.outofcontext.beta real number in the range [0;1] (default: 0.1).
#' Reads with average beta value for out-of-context cytosines \strong{above}
#' this threshold are considered completely unmethylated (all C are counted as
#' T). This option has no effect when read thresholding is disabled.
#' @param report.context string defining cytosine methylation context to report
#' (default: value of `threshold.context`).
#' @param ... other parameters to pass to the
#' \code{\link[epialleleR]{preprocessBam}} function.
#' Options have no effect if preprocessed BAM data was supplied as an input.
#' @param gzip boolean to compress the report (default: FALSE).
#' @param verbose boolean to report progress and timings (default: TRUE).
#' @return \code{\link[data.table]{data.table}} object containing cytosine
#' report in Bismark-like format or NULL if report.file was specified. The
#' report columns are:
#' \itemize{
#'   \item rname --- reference sequence name (as in BAM)
#'   \item strand --- strand
#'   \item pos --- cytosine position
#'   \item context --- methylation context
#'   \item meth --- number of methylated cytosines
#'   \item unmeth --- number of unmethylated cytosines
#' }
#' @seealso `values` vignette for a comparison and visualisation of epialleleR
#' output values for various input files. `epialleleR` vignette for the
#' description of usage and sample data.
#' 
#' \code{\link{preprocessBam}} for preloading BAM data,
#' \code{\link{generateBedReport}} for genomic region-based statistics,
#' \code{\link{generateVcfReport}} for evaluating epiallele-SNV associations,
#' \code{\link{extractPatterns}} for exploring methylation patterns,
#' \code{\link{generateBedEcdf}} for analysing the distribution of per-read
#' beta values.
#' @examples
#'   capture.bam <- system.file("extdata", "capture.bam", package="epialleleR")
#'   
#'   # CpG report with thresholding
#'   cg.report <- generateCytosineReport(capture.bam)
#'   
#'   # CX report without thresholding
#'   cx.report <- generateCytosineReport(capture.bam, threshold.reads=FALSE,
#'                report.context="CX")
#' @export
generateCytosineReport <- function (bam,
                                    report.file=NULL,
                                    threshold.reads=TRUE,
                                    threshold.context=c("CG", "CHG", "CHH", "CxG", "CX"),
                                    min.context.sites=2,
                                    min.context.beta=0.5,
                                    max.outofcontext.beta=0.1,
                                    report.context=threshold.context,
                                    ...,
                                    gzip=FALSE,
                                    verbose=TRUE)
{
  threshold.context <- match.arg(threshold.context, threshold.context)
  report.context    <- match.arg(report.context, report.context)
  
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
  
  cx.report <- .getCytosineReport(
    bam.processed=bam, pass=pass,
    ctx=.context.to.bases[[report.context]][["ctx.meth"]],
    verbose=verbose
  )
  
  if (is.null(report.file))
    return(cx.report)
  else
    .writeReport(report=cx.report, report.file=report.file, gzip=gzip,
                 verbose=verbose)
}

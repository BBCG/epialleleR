#' generateMhlReport
#'
#' @description
#' This function computes \emph{Linearised} Methylated Haplotype Load
#' (\eqn{lMHL}) per genomic position.
#'
#' @details
#' The function reports \emph{Linearised} Methylated Haplotype Load
#' (\eqn{lMHL}) at the level of individual cytosines using BAM file or data as
#' an input. Function uses the following formula:
#' 
#' \deqn{lMHL=\frac{\sum_{i=1}^{l} w_{i} \times MH_{i}}{\sum_{i=1}^{l} w_{i} \times H_{i}}}
#' 
#' where \eqn{l} is the length of a calculation window
#' (e.g., number of CpGs; \eqn{l \le L},
#' where \eqn{L} is the length of a haplotype covering current genomic
#' position),
#' \eqn{MH_{i}} is a number of fully successive methylated stretches
#' with \eqn{i} loci within fully successive methylated stretches that overlap
#' current genomic position,
#' \eqn{H_{i}} is a number of fully successive stretches with \eqn{i} loci,
#' \eqn{w_{i}} is a weight for \eqn{i}-locus haplotype (\eqn{w_{i}=i}).
#' 
#' This formula is a modification of the original 
#' Methylated Haplotype Load (MHL) formula that was 
#' first described by Guo et al., 2017
#' (doi: \href{https://doi.org/10.1038/ng.3805}{10.1038/ng.3805}): 
#' 
#' \deqn{MHL=\frac{\sum_{i=1}^{L} w_{i} \times P(MH_{i})}{\sum_{i=1}^{L} w_{i}}}
#' 
#' where \eqn{L} is the length of a longest haplotype covering current genomic
#' position, \eqn{P(MH_{i})=\frac{MH_{i}}{H_{i}}} is the fraction
#' of fully successive methylated stretches with \eqn{i} loci, \eqn{w_{i}} is
#' a weight for \eqn{i}-locus haplotype (\eqn{w_{i}=i}).
#' 
#' The modifications to original formula are made in order to:
#' \itemize{
#'   \item \bold{reduce the complexity of MHL calculation} for data of high
#'   breadth and depth --- \eqn{lMHL} values for all genomic positions can be
#'   calculated using a single pass (cycling through reads just once) as
#'   the linearised calculations of numerator and denominator for
#'   \eqn{lMHL} do not require prior knowledge on how many reads cover
#'   a particular position. This is achieved by moving \eqn{H_{i}}
#'   to the denominator of the \eqn{lMHL} formula.
#'   \item \bold{provide granularity of calculations/values?} --- the original MHL
#'   formula gives the same MHL value for every cytosine of a partially
#'   methylated haplotype (e.g., MHL=0.358 for each cytosine within a read with
#'   methylation call string "zZZZ"). 
#'   \item \bold{enable calculations for long-read seq}
#' }
#' 
#' These modifications do not affect ... but rather increase its value and
#' applicability for analysis of various sequencing data
#' 
#' However, this function introduces \bold{two optional parameters} that define
#' how \eqn{lMHL} formula is applied. These parameters are sought to make
#' \eqn{lMHL} calculations more flexible for various use case scenarios. These
#' parameters and scenarios are described below.
#' 
#' \bold{Case 1:}
#' 
#' 
#' 
#' 
#' 
#' In contrast to the other currently available software, reads
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
#' @param bam BAM file location string OR preprocessed output of
#' \code{\link[epialleleR]{preprocessBam}} function. BAM file alignment records
#' must contain XG tag (strand information for the reference genome) and
#' methylation call string (XM tag). Read more about these and other
#' requirements and BAM preprocessing at \code{\link{preprocessBam}}.
#' @param report.file file location string to write the cytosine report. If NULL
#' (the default) then report is returned as a
#' \code{\link[data.table]{data.table}} object.
#' @param haplotype.context string for a cytosine context that defines
#' a haplotype:
#' \itemize{
#'   \item "CG" (the default) -- CpG cytosines only (called as zZ)
#'   \item "CHG" -- CHG cytosines only (xX)
#'   \item "CHH" -- CHH cytosines only (hH)
#'   \item "CxG" -- CG and CHG cytosines (zZxX)
#'   \item "CX" -- all cytosines; this, as well as other non-CG contexts, may
#'   have little sense but still included for consistency
#' }
#' @param discontinuous boolean for discontinuous MHL calculation.
#' If FALSE (the default), classical MHL formula is used and MHL value for all
#' cytosines within read pair is the same (e.g., would equal 5/35 for each
#' cytosine of the read pair with methylation call string of "zZZzZ").
#' When TRUE, numerator of MHL equation varies according to the length of
#' fully methylated stretch the current cytosine belongs to (e.g., methylation
#' call string "zZZzZ" will produce the following MHL values per cytosine: 0/35, 
#' 4/35, 4/35, 0/35, 1/35). Discontinuous MHL calculations provide granularity
#' for reads / read pairs with span longer than expected haplotype blocks
#' (long-range sequencing) or for those that despite being short overlap with
#' more than one distinct haplotype block (amplicon sequencing). For thorough
#' explanation and more examples, see Details section and vignette.
#' @param max.l non-negative integer for maximum value of `l` in MHL formula.
#' When 0 (the default), MHL is calculated using classical formula where `l`
#' equals to the length of haplotype (number of CpG cytosines). When `max.l`>0,
#' `max.l` is used instead of `l` in MHL calculations. E.g., for two similar
#' read pairs
#' with methylation call strings of "ZZZZzZZZZ" and "ZZZZZZZZZ", classical MHL
#' values are 40/165 and 165/165, respectively. When calculated with `max.l`=4,
#' values are 40/70 and 70/70, respectively. Non-zero `max.l` values therefore
#' provide more granularity in MHL calculations and may help nivellate the
#' effect of sequencing errors for long reads. For thorough
#' explanation and more examples, see Details section and vignette.
#' @param min.L non-negative integer for minimum length of a haplotype
#' (default: 0). Reads (read pairs) with fewer than `min.L` cytosines
#' within the `haplotype.context` are skipped.
#' @param ... other parameters to pass to the
#' \code{\link[epialleleR]{preprocessBam}} function.
#' Options have no effect if preprocessed BAM data was supplied as an input.
#' @param gzip boolean to compress the report (default: FALSE).
#' @param verbose boolean to report progress and timings (default: TRUE).
#' @return \code{\link[data.table]{data.table}} object containing cytosine
#' report in Bismark-like format or NULL if report.file was specified. The
#' report columns are:
#' \itemize{
#'   \item rname -- reference sequence name (as in BAM)
#'   \item strand -- strand
#'   \item pos -- cytosine position
#'   \item context -- methylation context
#'   \item meth -- number of methylated cytosines
#'   \item unmeth -- number of unmethylated cytosines
#' }
#' @seealso \code{\link{preprocessBam}} for preloading BAM data,
#' \code{\link{generateBedReport}} for genomic region-based statistics,
#' \code{\link{generateVcfReport}} for evaluating epiallele-SNV associations,
#' \code{\link{extractPatterns}} for exploring methylation patterns,
#' \code{\link{generateBedEcdf}} for analysing the distribution of per-read
#' beta values, and `epialleleR` vignettes for the description of usage and
#' sample data.
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
generateMhlReport <- function (bam,
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

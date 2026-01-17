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
#' NB: you can modify and/or run this example -- see Examples section at
#' the bottom of this page.
#' 
#' Let's suppose we have a BAM file with four reads, all mapped to the "+"
#' strand of chromosome 1, positions 1-16. Assuming the default values
#' for the thresholding parameters (cytosine.context = "CG",
#' filter.reads=TRUE, max.outofcontext.beta = 0.1,
#' threshold.reads = TRUE, min.context.sites = 2, min.context.beta = 0.5),
#' the input and results will look as following:
#' 
#' \tabular{lllll}{
#'   methylation string \tab filter \tab threshold \tab explained \tab methylation reported \cr
#'   ...Z..x+.h..x..h. \tab excluded \tab <NA> \tab min.context.sites < 2 (only one zZ base) \tab all cytosines unmethylated \cr
#'   ...Z..z.h..x..h.  \tab pass \tab above \tab pass all criteria \tab only C4 (Z at position 4) is methylated \cr
#'   ...Z..z.h..X..h.  \tab excluded \tab <NA> \tab max.outofcontext.beta > 0.1 (1XH / 3xXhH = 0.33) \tab read excluded from reporting \cr
#'   ...Z..z.h..z-.h.  \tab pass \tab below \tab min.context.beta < 0.5 (1Z / 3zZ = 0.33) \tab all cytosines unmethylated
#' }
#' 
#' Since the reads number one and three are filtered out, and only the second
#' read will satisfy the thresholding criteria, the following CX report will be
#' produced (given that all reads map to chr1:+:1-16):
#' 
#' \tabular{llllll}{
#'   rname \tab strand \tab pos \tab context \tab meth \tab unmeth \cr
#'   chr1 \tab + \tab 4 \tab CG \tab 1 \tab 1 \cr
#'   chr1 \tab + \tab 7 \tab CG \tab 0 \tab 2 \cr
#'   chr1 \tab + \tab 9 \tab CHH \tab 0 \tab 2 \cr
#'   chr1 \tab + \tab 15 \tab CHH \tab 0 \tab 2 
#' }
#' 
#' (Disclaimer: the cytosine base at position 12 is absent in the report,
#' because its context cannot be determined based on two reads that passed the
#' filtering: in one of the reads this base is in CG context and in the other in
#' CHG context. Since none of these contexts is present in \strong{more} than a
#' half of the reads, the base is skipped.)
#' 
#' With the read filtering and thresholding disabled (filter.reads=FALSE,
#' threshold.reads = FALSE) all reads will be included and all methylated bases
#' will retain their status, so the CX report will be very similar
#' (nearly identical) to the reports produced by other methylation callers
#' (such as Bismark or Illumina DRAGEN Bio IT Platform):
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
#' @param cytosine.context string defining cytosine methylation context used
#' for filtering and/or thresholding the reads:
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
#' @param filter.reads boolean defining if sequence reads with too few context
#' bases or too high out-of-context cytosine methylation should be filtered
#' out (e.g., reads resulting from incompletely bisulfite-converted templates).
#' Default: TRUE.
#' @param min.context.sites non-negative integer for minimum number of cytosines
#' within the `cytosine.context` (default: 2). Reads containing \strong{fewer}
#' within-the-context cytosines will not be thresholded and will be ignored
#' in further computations.
#' This option has no effect when read filtering is disabled.
#' @param max.outofcontext.beta real number in the range [0;1] (default: 0.1).
#' Reads with average beta value for out-of-context cytosines \strong{above}
#' this threshold will not be thresholded and will be ignored in further
#' computations. This option has no effect when read filtering is disabled.
#' @param threshold.reads boolean defining if sequence reads (read pairs) should
#' be thresholded before counting methylated cytosines (default: TRUE).
#' Disabling thresholding (together with filtering)
#' makes the report virtually indistinguishable from the
#' ones generated by other software, such as Bismark or Illumina DRAGEN Bio IT
#' Platform. Thresholding is \strong{not} recommended for long-read sequencing
#' data because long reads might cover multiple regions with very different
#' cytosine methylation.
#' @param min.context.beta real number in the range [0;1] (default: 0.5). Reads
#' with average beta value for within-the-context cytosines \strong{below} this
#' threshold are considered completely unmethylated (all C are counted as T).
#' This option has no effect when read thresholding is disabled.
#' @param report.context string defining cytosine methylation context to report
#' (default: value of `cytosine.context`).
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
#' \code{\link{extractPatterns}} for exploring methylation patterns and
#' \code{\link{plotPatterns}} for pretty plotting of its output,
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
#'   
#'   # Long-read sequencing with both filtering and thresholding disabled
#'   long.bam <- system.file("extdata", "longread.bam", package="epialleleR")
#'   long.data <- preprocessBam(bam=long.bam, min.mapq=30, min.baseq=20,
#'                              min.prob=178)
#'   cg.report <- generateCytosineReport(bam=long.data, filter.reads=FALSE,
#'                                       threshold.reads=FALSE)
#'   plot(cg.report[, .(pos, beta=data.table::frollmean(meth/(meth+unmeth), 100))], type="l")
#'   
#'   # toy example from the description
#'   temp.bam <- tempfile(fileext=".bam") 
#'   temp.bed <- as("chr1:1-100", "GRanges")
#'   simulateBam(output.bam.file=temp.bam, rname="chr1", XG="CT",
#'               seq=c("AGACGTTAGTAATAGTA", "AAACGTTGTAATAGTA",
#'                     "AGACGTTGTAACAGTA",  "AAACGTTGTAATGTA"),
#'               XM=c( "...Z..x+.h..x..h.", "...Z..z.h..x..h.",
#'                     "...Z..z.h..X..h.",  "...Z..z.h..z.h."),
#'               cigar=c("7M1I9M", "16M", "16M", "12M1D3M"))
#'   # with read filtering and thresholding
#'   generateCytosineReport(bam=temp.bam, report.context="CX")
#'   # without read filtering
#'   generateCytosineReport(bam=temp.bam, report.context="CX",
#'                          filter.reads=FALSE)
#'   # without read thresholding
#'   generateCytosineReport(bam=temp.bam, report.context="CX",
#'                          threshold.reads=FALSE)
#'   # both read filtering and thresholding disabled = similar to other software
#'   generateCytosineReport(bam=temp.bam, report.context="CX",
#'                          filter.reads=FALSE, threshold.reads=FALSE)
#'   # patterns plotted
#'   plotPatterns(
#'     extractPatterns(bam=temp.bam, bed=temp.bed, extract.context="CX"),
#'     plot.context="CX", npatterns.per.bin=Inf
#'   )
#' @export
generateCytosineReport <- function (bam,
                                    report.file=NULL,
                                    cytosine.context=c("CG", "CHG", "CHH", "CxG", "CX"),
                                    filter.reads=TRUE,
                                    min.context.sites=2,
                                    max.outofcontext.beta=0.1,
                                    threshold.reads=TRUE,
                                    min.context.beta=0.5,
                                    report.context=cytosine.context,
                                    ...,
                                    gzip=FALSE,
                                    verbose=TRUE)
{
  cytosine.context <- match.arg(cytosine.context, cytosine.context)
  report.context   <- match.arg(report.context, report.context)
  
  bam <- preprocessBam(bam.file=bam, ..., verbose=verbose)
  
  pass <- .filterThresholdReads(
    bam.processed=bam,
    ctx.meth=.context.to.bases[[cytosine.context]][["ctx.meth"]],
    ctx.unmeth=.context.to.bases[[cytosine.context]][["ctx.unmeth"]],
    ooctx.meth=.context.to.bases[[cytosine.context]][["ooctx.meth"]],
    ooctx.unmeth=.context.to.bases[[cytosine.context]][["ooctx.unmeth"]],
    filter.reads=filter.reads,
    min.context.sites=min.context.sites,
    max.outofcontext.beta=max.outofcontext.beta,
    threshold.reads=threshold.reads,
    min.context.beta=min.context.beta,
    verbose=verbose
  )
  
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

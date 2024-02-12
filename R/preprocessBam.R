#' preprocessBam
#'
#' @description
#' This function reads and preprocesses BAM file.
#'
#' @details
#' The function loads and preprocesses BAM file, saving time when multiple
#' analyses are to be performed on large input files. Currently, HTSlib
#' is used to read the data, therefore it is possible to speed up the loading
#' by means of HTSlib decompression threads.
#' 
#' This function is also called internally when BAM file location is supplied as
#' an input for other `epialleleR` methods.
#' 
#' Please note that for BAM preprocessing as well as all its reporting methods,
#' `epialleleR` requires genomic
#' strand (XG tag) and a methylation call string (XM tag) to be present in a
#' BAM file - i.e., methylation calling must be performed after read
#' mapping/alignment by your software of choice. It is the case for BAM files
#' produced by Bismark Bisulfite Read Mapper and Methylation Caller,
#' Illumina DRAGEN, Illumina Cloud analysis solutions, as well as
#' contemporary Illumina sequencing instruments
#' with on-board read mapping/alignment (NextSeq 1000/2000, NovaSeq X),
#' therefore such files can be analysed without additional steps.
#' For alignments produced by other tools, e.g., BWA-meth, methylation calling
#' must be performed prior to BAM loading / reporting, by means of
#' \code{\link[epialleleR]{callMethylation}}.
#' 
#' `preprocessBam` always tests if BAM file is paired- or single-ended
#' and has all necessary tags (XM/XG) available. It is recommended to use
#' `verbose` processing and check messages for correct identification of
#' alignment endness. Otherwise, if the `paired` parameter is set explicitly,
#' exception is thrown when expected endness differs from the auto detected one.
#' 
#' During preprocessing of paired-end alignments, paired reads are merged
#' according to
#' their base quality: nucleotide base with the highest value in the QUAL string
#' is taken, unless its quality is less than `min.baseq`, which results in no
#' information for that particular position ("-"/"N"). These merged reads are
#' then processed as a single entity in all `epialleleR` methods. Due to
#' merging, overlapping bases in read pairs are counted only once, and the base
#' with the highest quality is taken.
#' 
#' During preprocessing of single-end alignments, no read merging is
#' performed. Only bases with quality of at least `min.baseq` are considered.
#' Lower base quality results in no information for that particular position
#' ("-"/"N").
#' 
#' For RRBS-like protocols, it is possible to trim alignments from one or both
#' ends. Trimming is performed during BAM loading and will therefore influence
#' results of all downstream `epialleleR` methods. Internally, trimming is
#' performed at the level of a template (i.e., read pair for paired-end BAM or
#' individual read for single-end BAM). This ensures that only necessary parts
#' (real ends of sequenced fragment) are removed for paired-end sequencing
#' reads.
#' 
#' It is also a requirement currently that paired-end BAM file must be sorted by
#' QNAME instead
#' of genomic location (i.e., "unsorted") to perform merging of paired-end
#' reads. Error message is shown if it is sorted by genomic location, in this
#' case please sort it by QNAME using 'samtools sort -n -o out.bam in.bam'.
#' 
#'
#' @param bam.file BAM file location string.
#' @param paired boolean for expected alignment endness: `TRUE` for paired-end,
#' `FALSE` for single-end, or `NULL` for auto detect (the default).
#' @param min.mapq non-negative integer threshold for minimum read mapping
#' quality (default: 0).
#' @param min.baseq non-negative integer threshold for minimum nucleotide base
#' quality (default: 0).
#' @param skip.duplicates boolean defining if duplicate aligned reads should be
#' skipped (default: FALSE). Option has no effect if duplicate reads were not
#' marked by alignment software.
#' @param trim non-negative integer or vector of length 2 for the number of
#' nucleotide bases to be trimmed from 5' and 3' ends of a template (i.e., 
#' read pair for paired-end BAM or read for single-end BAM).
#' Default: 0 for no trimming. Specifying `trim=1` will result in removing of
#' a single base from both ends, while specifying `trim=c(1,2)` will
#' result in removing of a single base from 5' end and 2 bases from 3' end.
#' @param nthreads non-negative integer for the number of additional HTSlib
#' threads to be used during BAM file decompression (default: 1). Two threads
#' (and usually no more than two) make sense for the files larger than 100 MB.
#' @param verbose boolean to report progress and timings (default: TRUE).
#' @return \code{\link[data.table]{data.table}} object containing preprocessed
#' BAM data.
#' @seealso \code{\link{preprocessGenome}} for preloading reference
#' sequences and \code{\link{callMethylation}} for methylation calling.
#' 
#' \code{\link{generateCytosineReport}} for methylation statistics at
#' the level of individual cytosines, \code{\link{generateBedReport}} for
#' genomic region-based statistics, \code{\link{generateVcfReport}} for
#' evaluating epiallele-SNV associations, \code{\link{extractPatterns}} for
#' exploring methylation patterns, \code{\link{generateBedEcdf}} for
#' analysing the distribution of per-read beta values, and `epialleleR`
#' vignettes for the description of usage and sample data.
#' 
#' Sequence Alignment/Map \href{https://samtools.github.io/hts-specs/SAMv1.pdf}{format specifications},
#' duplicate alignments marking by \href{http://www.htslib.org/doc/samtools-markdup.html}{Samtools}
#' and \href{https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/DuplicateMarking_fDG.htm}{Illumina DRAGEN Bio IT Platform}.
#' @examples
#'   capture.bam <- system.file("extdata", "capture.bam", package="epialleleR")
#'   bam.data    <- preprocessBam(capture.bam)
#' @export
preprocessBam <- function (bam.file,
                           paired=NULL,
                           min.mapq=0,
                           min.baseq=0,
                           skip.duplicates=FALSE,
                           trim=0,
                           nthreads=1,
                           verbose=TRUE)
{
  if (is.character(bam.file)) {
    paired.check <- .checkBam(bam.file=bam.file, verbose=verbose)
    if (!is.null(paired))
      if (paired.check!=paired)
        stop("Expected endness is different from detected! Exiting",call.=FALSE)
    trim <- utils::head(rep.int(trim, times=2), 2)
    bam.processed <- .readBam(
      bam.file=bam.file, paired=paired.check,
      min.mapq=min.mapq, min.baseq=min.baseq,
      skip.duplicates=skip.duplicates, trim=trim,
      nthreads=nthreads, verbose=verbose
    )
    return(bam.processed)
  } else {
    if (verbose &
        !all(missing(paired), missing(min.mapq), missing(min.baseq),
             missing(skip.duplicates), missing(trim), missing(nthreads))) 
      message("Already preprocessed BAM supplied as an input. Explicitly set",
              " 'preprocessBam' options will have no effect.")
    return(bam.file)
  }
}

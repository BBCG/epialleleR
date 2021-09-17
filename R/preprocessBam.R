#' preprocessBam
#'
#' @description
#' This function reads and preprocesses BAM file.
#'
#' @details
#' The function loads and preprocesses BAM file, saving time when multiple
#' analyses are to be performed on large input files. Currently, HTSlib
#' is used to read the data, therefore it is possible to speed up the loading
#' by means of HTSlib threads.
#' 
#' This function is also called internally when BAM file location is supplied as
#' an input for other `epialleleR` methods.
#' 
#' `preprocessBam` currently accepts only BAM files that are derived from
#' paired-end sequencing (create an issue if you need to process single-end BAM
#' files). During preprocessing, paired reads are merged according to their base 
#' quality: nucleotide base with the highest value in the QUAL string is taken,
#' unless its quality is less than min.baseq, which results in no information
#' for that particular position ("-"/"N"). These merged reads are then
#' processed as a single entity in all `epialleleR` methods.
#' 
#' Please also note that for all its methods, `epialleleR` requires methylation
#' call string to be present in a BAM file - i.e., methylation calling must be
#' performed after read mapping/alignment by your software of choice.
#'
#' @param bam.file BAM file location string.
#' @param min.mapq non-negative integer threshold for minimum read mapping
#' quality (default: 0).
#' @param min.baseq non-negative integer threshold for minimum nucleotide base
#' quality (default: 0).
#' @param skip.duplicates boolean defining if duplicate aligned reads should be
#' skipped (default: FALSE). Option has no effect if duplicate reads were not
#' marked by alignment software.
#' @param nthreads non-negative integer for the number of HTSlib threads to be
#' used during BAM file decompression (default: 1). 2 threads make sense for the
#' files larger than 100 MB.
#' @param verbose boolean to report progress and timings (default: TRUE).
#' @return \code{\link[data.table]{data.table}} object containing preprocessed
#' BAM data.
#' @seealso \code{\link{generateCytosineReport}} for methylation statistics at
#' the level of individual cytosines, \code{\link{generateBedReport}} for
#' genomic region-based statistics, \code{\link{generateVcfReport}} for
#' evaluating epiallele-SNV associations, \code{\link{generateBedEcdf}} for
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
                           min.mapq=0,
                           min.baseq=0,
                           skip.duplicates=FALSE,
                           nthreads=1,
                           verbose=TRUE)
{
  if (is.character(bam.file)) {
    bam.processed <- .readBam(
      bam.file=bam.file, min.mapq=min.mapq, min.baseq=min.baseq,
      skip.duplicates=skip.duplicates, nthreads=nthreads, verbose=verbose
    )
    return(bam.processed)
  } else {
    if (verbose) 
      message("Already preprocessed BAM supplied as an input. Options",
              " 'min.mapq', 'min.baseq', 'skip.duplicates' and 'nthreads' ",
              "will have no effect.")
    return(bam.file)
  }
}

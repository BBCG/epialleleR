#' preprocessBam
#'
#' @description
#' This function reads and preprocesses BAM file.
#'
#' @details
#' The function loads and preprocesses BAM file, saving time when multiple
#' analyses are to be performed on large input files. Currently, Rsamtools
#' package is used to read the data, but this will change in the future with a
#' goal of speeding up this step even further.
#' 
#' This function is also used internally when BAM file location is supplied as
#' an input for other `epialleleR` methods.
#' 
#' `preprocessBam` automatically determines whether BAM is derived from
#' single-end or paired-end sequencing. When the latter is the case, paired
#' reads are merged so that the overlapping fragments of second read are clipped
#' (because quality of the second read is usually lower than of the first).
#' These merged reads are then processed as a single entity in all `epialleleR`
#' methods.
#' 
#' Please also note that for all its methods, `epialleleR` requires methylation
#' call string to be present in a BAM file - i.e. methylation calling must be
#' performed after read mapping/alignment by your software of choice.
#'
#' @param bam.file BAM file location string.
#' @param min.mapq non-negative integer threshold for minimum read mapping
#' quality (default: 0).
#' @param skip.duplicates boolean defining if duplicate aligned reads should be
#' skipped (default: FALSE). Option has no effect if duplicate reads were not
#' marked by alignment software.
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
                           skip.duplicates=FALSE,
                           verbose=TRUE)
{
  if (is.character(bam.file)) {
    bam <- .readBam(bam.file=bam.file, min.mapq=min.mapq,
                    skip.duplicates=skip.duplicates, verbose=verbose)
    bam.processed <- .processBam(bam=bam, verbose=verbose)
    return(bam.processed)
  } else {
    if (verbose) 
      message("Already preprocessed BAM supplied as an input. Options",
              " 'min.mapq' and 'skip.duplicates' will have no effect.")
    return(bam.file)
  }
}

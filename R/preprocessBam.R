#' preprocessBam
#'
#' @description
#' `preprocessBam` This function reads and preprocesses BAM file
#'
#' @details
#' The function loads and preprocesses BAM file, saving time when multiple
#' analyses are to be performed on large input files. Currently, Rsamtools
#' package is used to read the data, but this will change in a future with a
#' goal of speeding up this step even further.
#' 
#' This function is also used internally when BAM file location is supplied as
#' an input for other `epialleleR` methods
#' 
#' Please also note that for all its methods, `epialleleR` requires methylation
#' call string to be present in a BAM file - i.e. methylation calling must be
#' performed after read mapping/alignment by your software of choice.
#'
#' @param bam.file BAM file location string
#' @param min.mapq non-negative integer threshold for minimum read mapping
#' quality (default: 0)
#' @param skip.duplicates boolean defining if duplicate aligned reads should be
#' skipped (default: FALSE). Option has no effect if duplicate reads were not
#' marked by alignment software
#' @param verbose boolean to report progress and timings (default: TRUE)
#' @return \code{\link[data.table]{data.table}} object containing preprocessed
#' BAM data
#' @seealso \href{https://samtools.github.io/hts-specs/SAMv1.pdf}{Sequence Alignment/Map format specifications},
#' duplicate alighments marking by \href{http://www.htslib.org/doc/samtools-markdup.html}{Samtools}
#' and \href{https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/DuplicateMarking_fDG.htm}{Illumina DRAGEN Bio IT Platform}
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


#|[c]{}^

# ##############################################################################

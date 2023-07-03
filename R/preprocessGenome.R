#' preprocessGenome
#'
#' @description
#' This function reads and preprocesses (optionally `bgzip`ped) FASTA file with
#' reference sequences.
#'
#' @details
#' The function loads and preprocesses reference (genomic) sequences, saving
#' time when methylation calling needs to be performed on multiple BAM files.
#' Currently, reading the data is done by means of HTSlib,
#' therefore it is possible to speed up the loading
#' by means of HTSlib decompression threads when FASTA file is compressed by
#' `bgzip`. 
#' 
#' This function is also called internally when file location is supplied
#' as an input for \code{\link{callMethylation}} method.
#' 
#' `preprocessGenome` checks if index file is present, and if not, creates
#' it automatically. It is possible and recommended to use compressed FASTA
#' file as an input, but the file must be compressed by `bgzip` (part of
#' samtools/HTSlib). When FASTA file is compressed, faster loading can be
#' achieved using (typically one) additional HTSlib decompression thread.
#' 
#' Please also note that for the purpose of methylation calling, the very same
#' reference genome must be used for both alignment (when BAM is produced) and
#' calling cytosine methylation by \code{\link{callMethylation}} method.
#'
#' @param genome.file reference (genomic) sequences file location string.
#' @param nthreads non-negative integer for the number of additional HTSlib
#' threads to be used during file decompression (default: 1).
#' @param verbose boolean to report progress and timings (default: TRUE).
#' @return list object containing preprocessed reference sequence data.
#' @seealso \code{\link{callMethylation}} for methylation calling,
#' and `epialleleR` vignettes for the description of usage and sample data.
#' 
#' Block compression/decompression utility \href{http://www.htslib.org/doc/bgzip.html}{bgzip}.
#' @examples
#'   genome.file <- system.file("extdata", "test", "reference.fasta.gz", package="epialleleR")
#'   genome.data <- preprocessGenome(genome.file)
#' @export
preprocessGenome <- function (genome.file,
                              nthreads=1,
                              verbose=TRUE)
{
  if (is.character(genome.file)) {
    genome.processed <- .readGenome(
      genome.file=genome.file, nthreads=nthreads, verbose=verbose
    )
    return(genome.processed)
  } else {
    return(genome.file)
  }
}

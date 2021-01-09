#' preprocessBam
#'
#' @description
#' `preprocessBam` desc. Epiallele-aware
#'
#' @details
#' details
#'
#' @param bam.file param desc
#' @param min.mapq param desc
#' @param skip.duplicates param desc
#    http://www.htslib.org/doc/samtools-markdup.html
#    https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/DuplicateMarking_fDG.htm
#' @param verbose param desc
#' @return return desc
#' @seealso whatelse
#' @examples
#' \dontrun{
#'   sessionInfo()
#' }
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
  }

  return(bam.processed)
}


#|[c]{}^

# ##############################################################################

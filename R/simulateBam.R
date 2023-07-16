#' simulateBam
#'
#' @description
#' This function creates sample BAM files given mandatory and optional BAM
#' fields.
#'
#' @details
#' The function creates sample alignment records and saves them in BAM file.
#' Output can be used to test epialleleR methods as well as other
#' tools for methylation analysis. This method can significantly simplify
#' calculation of methylation metrics on example data (beta, VEF, and WTF
#' values of epialleleR; methylation heterogeneity metrics of other tools).
#' 
#' The number of records written will be equal to the largest length of any
#' supplied (nondefault) parameter or 1 if no parameters were supplied.
#' If lengths of supplied parameters differ,
#' shorter vectors will be recycled (a whole number of times or with remainder
#' if necessary).
#' 
#' Please note that function performs almost no validity checks for supplied
#' fields. In particular, be extra careful constructing paired-end BAM
#' alignments, and if necessary use `samtools` to perform validity check or
#' manual editing after BAM->SAM conversion.
#'
#' @param output.bam.file output BAM file location string. If NULL (default),
#' records are not written to BAM but returned as a
#' \code{\link[data.table]{data.table}} object for review.
#' @param qname character vector of query names. When default (NULL), names
#' like "q0001".."qnnnn" will be assigned.
#' @param flag integer vector of bitwise flags (a combination of the BAM_F*
#' constants). When default (NULL), zero (i.e., unique, valid, single-end,
#' aligned read) is assigned for every record.
#' @param rname character vector of chromosome (reference) names. When default
#' (NULL), "chrS" is assigned for every record.
#' @param pos integer vector of 1-based leftmost coordinates of the queries.
#' When default (NULL), 1 is assigned for every record.
#' @param mapq integer vector of mapping qualities. When default (NULL),
#' 60 is assigned for every record.
#' @param cigar character vector of CIGAR strings. When default (NULL),
#' "lM" is assigned for every record, where `l` is the length of the query
#' (`seq`).
#' @param rnext character vector of chromosome (reference) names for next read
#' in template. When default (NULL), "chrS" is assigned for every record.
#' @param pnext integer vector of 1-based leftmost coordinates of next read in
#' template. When default (NULL), 1 is assigned for every record.
#' @param tlen integer vector of observed template lengths. When default
#' (NULL), the length of the corresponding query (`seq`)
#' is assigned for every record.
#' @param seq character vector of query sequences. When default (NULL),
#' random sequence is assigned.
#' The lengths of these random sequences equal to the lengths of
#' methylation call strings from the `XM` optional parameter (if supplied),
#' or to the `tlen` parameter (if defined).
#' If none of these parameters is supplied, length of every `seq` will equal 10.
#' @param qual query sequence quality strings (ASCII of base QUALity plus 33).
#' When default (NULL), quality of every base is assigned to "F" (QUALity
#' of 47 + 33). The lengths of these quality strings equal to the length of the
#' corresponding query sequences (`seq`) for every record.
#' @param ... optional tags to add to the records, in the form `tag=value`.
#' Can be either integer vector (e.g., for "NM" tag),
#' or character vector (e.g., "XM" tag for methylation call string,
#' "XG"/"YC" tag for reference strand read was aligned to).
#' @param verbose boolean to report progress and timings (default: TRUE).
#' @return number of BAM records written (if `output.bam.file` is not NULL) or
#' \code{\link[data.table]{data.table}} object containing final records
#' prepared for writing. NB: this object has 0-based coordinates and
#' numerically encoded reference names.
#' @seealso \code{\link{generateCytosineReport}} and
#' \code{\link{generateWtfReport}} for methylation reports, as well as
#' `epialleleR` vignettes for the description of usage and sample data.
#' 
#' \href{https://www.htslib.org/}{Samtools} for viewing BAM files.
#' \href{http://samtools.github.io/hts-specs/SAMv1.pdf}{SAMv1} file format
#' specifications. Specifications of
#' \href{https://samtools.github.io/hts-specs/SAMtags.pdf}{optional SAM tags}.
#' \href{https://doi.org/10.1371/journal.pcbi.1010946}{metheor} for ultrafast
#' DNA methylation heterogeneity calculation from bisulfite alignments.
#' @examples
#'   out.bam <- tempfile(pattern="simulated", fileext=".bam")
#'   simulateBam(
#'     output.bam.file=out.bam,
#'     pos=c(1, 2),
#'     XM=c("ZZZzzZZZ", "ZZzzzzZZ"),
#'     XG=c("CT", "AG")
#'   )
#'   generateCytosineReport(out.bam, threshold.reads=FALSE)
#' @export
simulateBam <- function (output.bam.file=NULL,
                         qname=NULL,
                         flag=NULL,
                         rname=NULL,
                         pos=NULL,
                         mapq=NULL,
                         cigar=NULL,
                         rnext=NULL,
                         pnext=NULL,
                         tlen=NULL,
                         seq=NULL,
                         qual=NULL,
                         ...,
                         verbose=TRUE)
{
  result <- .simulateBam(output.bam.file=output.bam.file,
                         qname=qname, flag=flag, rname=rname,
                         pos=pos, mapq=mapq, cigar=cigar, rnext=rnext,
                         pnext=pnext, tlen=tlen, seq=seq, qual=qual, ...,
                         verbose=verbose)
  return(result)
}

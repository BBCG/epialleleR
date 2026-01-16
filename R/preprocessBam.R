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
#' `preprocessBam` currently allows to load both short-read (e.g., bisulfite)
#' and long-read (native) sequencing alignments. Specific requirements for these
#' types of data are given below.
#' 
#' @section Short-read sequencing:
#' 
#' For preprocessing of short reads (and therefore for all reporting
#' methods), `epialleleR` requires genomic
#' strand (XG tag) and a methylation call string (XM tag) to be present in a
#' BAM file - i.e., methylation calling must be performed after read
#' mapping/alignment by your software of choice. It is the case for BAM files
#' produced by Bismark Bisulfite Read Mapper and Methylation Caller,
#' Illumina DRAGEN, Illumina Cloud analysis solutions, as well as
#' contemporary Illumina sequencing instruments
#' with on-board read mapping/alignment (NextSeq 1000/2000, NovaSeq X),
#' therefore such files can be analysed without additional steps.
#' For alignments produced by other tools, e.g., BWA-meth or BSMAP, methylation
#' calling must be performed prior to BAM loading / reporting, by means of
#' \code{\link[epialleleR]{callMethylation}}.
#' 
#' @section Long-read sequencing:
#' 
#' For preprocessing of long reads, `epialleleR` requires presence of MM (Mm)
#' and ML (Ml) tags that hold information on base modifications and related
#' probabilities, respectively. These are standard tags described in SAM/BAM
#' format specification, therefore relevant tools for analysis and alignment
#' of long sequencing reads should be able to produce them. 
#' 
#' @section Other details:
#' 
#' `preprocessBam` always tests if BAM file is paired- or single-ended
#' and has all the necessary tags available. It is recommended to use
#' `verbose` processing and check messages for correct identification of
#' alignment endness. Otherwise, if the `paired` parameter is set explicitly,
#' exception is thrown when expected endness differs from the autodetected one.
#' It is nevertheless possible to override the autodetected endness and load
#' BAM as specified in `paired` by setting the `override.check` to TRUE.
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
#' It is currently a requirement that paired-end BAM file must be sorted by
#' QNAME instead
#' of genomic location (i.e., "unsorted") to perform merging of paired-end
#' reads. Error message is shown if it is sorted by genomic location, in this
#' case please sort it by QNAME using 'samtools sort -n -o out.bam in.bam'.
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
#' It is also possible to load only a subset of reads (read pairs) of interest
#' or only fragments of such reads by supplying
#' a list of targets (see description of `targets` and other related options).
#' The subsetting is performed without using BAM index.
#' 
#' @section Specific considerations for long-read sequencing data:
#' 
#' Any location not reported is implicitly assumed to contain no modification.
#' 
#' According to SAM format specification,
#' MM base modification tags are allowed to list modifications observed not
#' only on the original sequenced strand (e.g., `C+m`) but also on the 
#' opposite strand (e.g., `G-m`). The logic of their processing is as follows
#' (with the examples given below):
#' 
#' \itemize{
#'   \item if an alignment record has no methylation modifications (neither
#'   `C+m`, nor `G-m` are present), this record is, naturally, considered to
#'   be a single read with no cytosines methylated
#'   \item if an alignment record has `C+m` modification (base modifications
#'   on the original sequenced strand), then this record is, naturally,
#'   considered to be a single read with cytosine modifications on the
#'   sequenced strand
#'   \item if an alignment record has `G-m` modification (base modifications
#'   on the strand opposite to sequenced), then this record is treated as two
#'   reads, with the original sequenced strand having no modifications,
#'   while the opposite strand having cytosine modifications
#'   \item if both `C+m` and `G-m` are present, then this record is treated
#'   as two reads, with both strands having cytosine modifications
#' }
#' 
#' 
#' @param bam.file BAM file location string.
#' @param paired boolean for expected alignment endness: `TRUE` for paired-end,
#' `FALSE` for single-end, or `NULL` for auto detect (the default).
#' @param override.check boolean to use supplied endness (`paired` parameter)
#' even if it is different from the autodetected one (default: FALSE).
#' @param min.mapq non-negative integer threshold for minimum read mapping
#' quality. Default is 0, however a higher value (e.g., 30) is recommended.
#' @param min.baseq non-negative integer threshold for minimum nucleotide base
#' quality. Default is 0, however a higher value (e.g., at least 13 as in
#' `samtools mpileup`) is recommended.
#' @param min.prob integer threshold for minimum scaled probability of
#' modification (methylation) to consider. Affects processing of long-read
#' sequencing alignments only. According to SAM/BAM specification, the
#' continuous base modification probability range 0.0 to 1.0 is
#' remapped in equal sized portions to the discrete integers 0 to 255
#' inclusively. Default is -1, however a higher value (e.g., 178 which
#' corresponds to a modification probability of ~0.7) is recommended.
#' Also, when default (-1), all C+m and G-m cytosine
#' methylation modifications recorded in MM/Mm tag will be included, even if
#' ML/Ml tag with probabilities is absent (in such case, probability of
#' modification equals -1).
#' @param highest.prob boolean defining if methylation modification must have
#' the highest probability among all modifications at a particular base to be
#' considered in the analyses (default: TRUE). Affects processing long-read
#' sequencing alignments only. If default (TRUE) and ML/Ml tag with probability
#' scores is absent, then cytosines with more than one modification will be
#' omitted (as the probability of all modifications will be equal).
#' @param skip.duplicates boolean defining if duplicate aligned reads should be
#' skipped (default: FALSE). Option has no effect if duplicate reads were not
#' marked by alignment software.
#' @param skip.secondary boolean defining if secondary alignments should be
#' skipped (default: TRUE). Do not change.
#' @param skip.qcfail boolean defining if alignments failing QC should be
#' skipped (default: TRUE). Do not change.
#' @param skip.supplementary boolean defining if supplementary alignments
#' should be skipped (default: TRUE). Do not change.
#' @param trim non-negative integer or vector of length 2 for the number of
#' nucleotide bases to be trimmed from 5' and 3' ends of a template (i.e., 
#' read pair for paired-end BAM or read for single-end BAM).
#' Default: 0 for no trimming. Specifying `trim=1` will result in removing of
#' a single base from both ends, while specifying `trim=c(1,2)` will
#' result in removing of a single base from 5' end and 2 bases from 3' end.
#' @param targets Browser Extensible Data (BED) file location string OR object
#' of class \code{\link[GenomicRanges]{GRanges}} holding genomic coordinates for
#' regions of interest. The style of seqlevels of BED file/object must be the
#' same as the style of seqlevels of BAM file. During targets loading, the 
#' genomic regions are reduced (i.e., overlapping regions are merged),
#' and the strand information is discarded. These reduced regions are then used
#' to only load overlapping reads (merged read pairs for paired-end sequencing).
#' @param zero.based.targets boolean defining if BED file coordinates are
#' zero-based (default: FALSE).
#' @param clip.to.targets boolean defining if overlapping reads (merged read
#' pairs in case of paired-end sequencing) should be clipped to retain only
#' fragments overlapping with `targets` (default: FALSE). If a read overlaps
#' with two or more targets, all overlapping fragments will be retained.
#' @param nthreads non-negative integer for the number of additional HTSlib
#' threads to be used during BAM file decompression (default: 1). Two threads
#' (and usually no more than four) make sense for files larger than 100 MB.
#' @param verbose boolean to report progress and timings (default: TRUE).
#' @return \code{\link[data.table]{data.table}} object containing preprocessed
#' BAM data.
#' 
#' NB: most of the BAM data is stored not in the data.table
#' \emph{per se}, but as an external object linked to this data.table.
#' Therefore, saving/loading this data.table cannot
#' be used to save/recover preprocessed BAM data.
#' @seealso \code{\link{preprocessGenome}} for preloading reference
#' sequences and \code{\link{callMethylation}} for methylation calling.
#' 
#' \code{\link{generateCytosineReport}} for methylation statistics at
#' the level of individual cytosines, \code{\link{generateBedReport}} for
#' genomic region-based statistics, \code{\link{generateVcfReport}} for
#' evaluating epiallele-SNV associations, \code{\link{extractPatterns}} for
#' exploring methylation patterns and \code{\link{plotPatterns}} for pretty
#' plotting of its output, \code{\link{generateBedEcdf}} for
#' analysing the distribution of per-read beta values, and `epialleleR`
#' vignettes for the description of usage and sample data.
#' 
#' Sequence Alignment/Map \href{https://samtools.github.io/hts-specs/SAMv1.pdf}{format specifications},
#' specifications for \href{https://samtools.github.io/hts-specs/SAMtags.pdf}{optional SAM tags},
#' duplicate alignments marking by \href{http://www.htslib.org/doc/samtools-markdup.html}{Samtools}
#' and \href{https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/DuplicateMarking_fDG.htm}{Illumina DRAGEN Bio IT Platform}.
#' @examples
#'   capture.bam <- system.file("extdata", "capture.bam", package="epialleleR")
#'   bam.data    <- preprocessBam(capture.bam)
#'   
#'   # Specifics of long-read alignment processing
#'   out.bam <- tempfile(pattern="out-", fileext=".bam")
#'   
#'   simulateBam(
#'     seq=c("ACGCCATYCGCGCCA"),
#'     Mm=c("C+m,0,2,0;"),
#'     Ml=list(as.integer(c(102,128,153))),
#'     output.bam.file=out.bam
#'   )
#'   generateCytosineReport(out.bam, threshold.reads=FALSE, report.context="CX")
#'   
#'   simulateBam(
#'     seq=c("ACGCCATYCGCGCCA"),
#'     Mm=c("G-m,0,0,0;"),
#'     Ml=list(as.integer(c(138,101,96))),
#'     output.bam.file=out.bam
#'   )
#'   generateCytosineReport(out.bam, threshold.reads=FALSE, report.context="CX")
#'   
#'   simulateBam(
#'     seq=c("ACGCCATYCGCGCCA"),
#'     Mm=c("C+m,0,2,0;G-m,0,0,0;"),
#'     Ml=list(as.integer(c(102,128,153,138,101,96))),
#'     output.bam.file=out.bam
#'   )
#'   generateCytosineReport(out.bam, threshold.reads=FALSE, report.context="CX")
#'   
#' @export
preprocessBam <- function (bam.file,
                           paired=NULL,
                           override.check=FALSE,
                           min.mapq=0,
                           min.baseq=0,
                           min.prob=-1,
                           highest.prob=TRUE,
                           skip.duplicates=FALSE,
                           skip.secondary=TRUE,
                           skip.qcfail=TRUE,
                           skip.supplementary=TRUE,
                           trim=0,
                           targets=NULL,
                           zero.based.targets=FALSE,
                           clip.to.targets=FALSE,
                           nthreads=1,
                           verbose=TRUE)
{
  if (is.character(bam.file)) {
    bam.check <- .checkBam(bam.file=bam.file, verbose=verbose)
    if (!is.null(paired) && bam.check$paired!=paired) {
      if (override.check) {
        warning("Supplied endness is different from detected!", call.=FALSE)
        bam.check$paired <- paired
      } else {
        stop("Supplied endness is different from detected! If it is a mistake,",
             " override with 'override.check=TRUE'", call.=FALSE)
      }
    }
    
    if (!is.null(targets) && !methods::is(targets, "GRanges"))
      targets <- .readBed(bed.file=targets, zero.based.bed=zero.based.targets,
                          verbose=verbose)
      
    trim <- rep_len(trim, length.out=2)
    bam.processed <- .readBam(
      bam.file=bam.file, bam.check=bam.check,
      min.mapq=min.mapq, min.baseq=min.baseq,
      min.prob=min.prob, highest.prob=highest.prob,
      skip.duplicates=skip.duplicates, skip.secondary=skip.secondary,
      skip.qcfail=skip.qcfail, skip.supplementary=skip.supplementary,
      trim=trim, targets=targets, clip.to.targets=clip.to.targets,
      nthreads=nthreads, verbose=verbose
    )
    return(bam.processed)
  } else {
    if (verbose &&
        !all(missing(paired), missing(override.check),
             missing(min.mapq), missing(min.baseq),
             missing(min.prob), missing(highest.prob),
             missing(skip.duplicates), missing(skip.secondary),
             missing(skip.qcfail), missing(skip.supplementary),
             missing(trim), missing(targets), missing(zero.based.targets),
             missing(clip.to.targets), missing(nthreads))) 
      message("Already preprocessed BAM supplied as an input. Explicitly set",
              " 'preprocessBam' options will have no effect.")
    return(bam.file)
  }
}

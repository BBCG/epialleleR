#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table as.data.table
#' @importFrom data.table dcast
#' @importFrom data.table merge.data.table
#' @importFrom data.table setorder
#' @importFrom data.table setkey
#' @importFrom data.table setDT
#' @importFrom data.table setattr
#' @importFrom stringi stri_length
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges reduce
#' @importFrom BiocGenerics sort
#' @importFrom BiocGenerics width
#' @importFrom VariantAnnotation ScanVcfParam
#' @importFrom VariantAnnotation readVcf
#' @importFrom VariantAnnotation expand
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom stats ecdf
#' @importFrom stats setNames
#' @importFrom methods is
#' @importFrom utils globalVariables
#' @useDynLib epialleleR, .registration=TRUE


# internal globals, constants and helper functions 
#

################################################################################
# Globals, unload
################################################################################

utils::globalVariables(
  c(".", ".I", ".N", ":=", "FALSE+", "FALSE-", "TRUE+", "TRUE-", "XG", "XM",
    "XM.norm", "bedmatch", "cigar", "context", "flag", "isize", "meth", "mpos",
    "pass", "pos", "qname", "rname", "strand", "triad", "unmeth", "width",
    "isfirst", "ALT", "M+A", "M+Alt", "M+C", "M+G", "M+Ref", "M+T", "M-A",
    "M-Alt", "M-C", "M-G", "M-Ref", "M-T", "REF", "U+A", "U+Alt", "U+C", "U+G",
    "U+Ref", "U+T", "U-A", "U-Alt", "U-C", "U-G", "U-Ref", "U-T")
)

.onUnload <- function (libpath) {library.dynam.unload("epialleleR", libpath)}

################################################################################
# Constants
################################################################################

# descr: methylation context to XM bases.
#        The "u" and "U" are simply ignored - just as Bismark does
#        in order to save us from confusion at the analysis stage

.context.to.bases <- list (
  "CG"  = list (ctx.meth   = "Z",   ctx.unmeth   = "z",
                ooctx.meth = "XH",  ooctx.unmeth = "xh"),
  "CHG" = list (ctx.meth   = "X",   ctx.unmeth   = "x",
                ooctx.meth = "ZH",  ooctx.unmeth = "zh"),
  "CHH" = list (ctx.meth   = "H",   ctx.unmeth   = "h",
                ooctx.meth = "ZX",  ooctx.unmeth = "zx"),
  "CxG" = list (ctx.meth   = "ZX",  ctx.unmeth   = "zx",
                ooctx.meth = "H",   ooctx.unmeth = "h"),
  "CX"  = list (ctx.meth   = "ZXH", ctx.unmeth   = "zxh",
                ooctx.meth = "",    ooctx.unmeth = "")
)

################################################################################
# Functions: reading/writing files
################################################################################

# descr: reads BAM files using HTSlib
# value: data.table

.readBam <- function (bam.file,
                      min.mapq,
                      min.baseq,
                      skip.duplicates,
                      nthreads,
                      verbose)
{
  if (verbose) message("Reading BAM file", appendLF=FALSE)
  tm <- proc.time()
  
  bam.file <- path.expand(bam.file)
  bam.processed <- rcpp_read_bam_paired(bam.file, min.mapq, min.baseq, 
                                        skip.duplicates, nthreads)
  data.table::setDT(bam.processed)
  data.table::setorder(bam.processed, rname, start)
  
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(bam.processed)
}

################################################################################

# descr: (fast) reads BED file with amplicons
# value: object of type GRanges

.readBed <- function (bed.file,
                      zero.based.bed,
                      verbose)
{
  if (verbose) message("Reading BED file", appendLF=FALSE)
  tm <- proc.time()
  
  bed.df <- data.table::fread(file=bed.file, sep="\t", blank.lines.skip=TRUE,
                              data.table=FALSE)
  colnames(bed.df)[seq_len(3)] <- c("chr", "start", "end")
  bed <- GenomicRanges::makeGRangesFromDataFrame(
    bed.df, keep.extra.columns=TRUE, ignore.strand=TRUE,
    starts.in.df.are.0based=zero.based.bed
  )
  
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(bed)
}

################################################################################

# descr: (fast) reads VCF file
# value: object of VCF type (VariantAnnotation)

.readVcf <- function (vcf.file,
                      vcf.style,
                      bed,
                      verbose)
{
  if (verbose) message("Reading VCF file", appendLF=FALSE)
  tm <- proc.time()
  
  vcf.genome <- "unknown"
  vcf.param  <- VariantAnnotation::ScanVcfParam(info=NA, geno=NA)
  
  if (!is.null(bed)) {
    bed.style  <- GenomeInfoDb::seqlevelsStyle(bed)
    if (!is.null(vcf.style))
      GenomeInfoDb::seqlevelsStyle(bed) <- vcf.style
    bed.levels <- levels(GenomicRanges::seqnames(bed))
    vcf.genome <- stats::setNames(rep("unknown",length(bed.levels)), bed.levels)
    vcf.param  <- VariantAnnotation::ScanVcfParam(info=NA, geno=NA, which=bed)
  }
  
  vcf <- VariantAnnotation::readVcf(
    file=vcf.file,
    genome=vcf.genome,
    param=vcf.param
  )
  
  if (exists("bed.style")) {
    vcf.ranges <- SummarizedExperiment::rowRanges(vcf)
    suppressPackageStartupMessages({
      GenomeInfoDb::seqlevelsStyle(vcf.ranges) <- bed.style
    })
    SummarizedExperiment::rowRanges(vcf) <- vcf.ranges
  }
  
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(vcf)
}

################################################################################

# descr: (fast) writes the report
# value: void

.writeReport <- function (report,
                          report.file,
                          gzip,
                          verbose)
{
  if (verbose) message("Writing the report", appendLF=FALSE)
  tm <- proc.time()
  
  data.table::fwrite(report, file=report.file, quote=FALSE, sep="\t",
                     row.names=FALSE, col.names=TRUE,
                     compress=if (gzip) "gzip" else "none")
  
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
}

################################################################################
# Functions: processing
################################################################################

# descr: apply thresholding criteria to processed BAM reads
# value: bool vector with true for reads passing the threshold

.thresholdReads <- function (bam.processed,
                             ctx.meth, ctx.unmeth, ooctx.meth, ooctx.unmeth,
                             min.context.sites, min.context.beta,
                             max.outofcontext.beta, verbose)
{
  if (verbose) message("Thresholding reads", appendLF=FALSE)
  tm <- proc.time()
  
  # fast thresholding, vectorised
  pass <- rcpp_threshold_reads(
    bam.processed,
    ctx.meth, ctx.unmeth, ooctx.meth, ooctx.unmeth,
    min.context.sites, min.context.beta, max.outofcontext.beta
  )
  
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(pass)
}

################################################################################

# descr: matching BED target (amplicon/capture)
# value: numeric vector

.matchTarget <- function (bam.processed, bed, bed.type,
                          match.tolerance, match.min.overlap)
{
  # fast, vectorised
  bed.dt <- data.table::as.data.table(bed)
  bed.dt[, seqnames := factor(seqnames, levels=levels(bam.processed$rname))]
  
  if (bed.type=="amplicon") {
    bed.match <- rcpp_match_amplicon(bam.processed, bed.dt, match.tolerance)
  } else if (bed.type=="capture") {
    bed.match <- rcpp_match_capture(bam.processed, bed.dt, match.min.overlap)
  }
  
  return(bed.match)
}

################################################################################
# Functions: reporting
################################################################################

# descr: prepare cytosine report for processed reads according to filter
# value: data.table with Bismark-like cytosine report

.getCytosineReport <- function (bam.processed,
                                pass,
                                ctx,
                                verbose)
{
  if (verbose) message("Preparing cytosine report", appendLF=FALSE)
  tm <- proc.time()
  
  # check if ordered? reorder if not
  # reporting, vectorised
  cx.report <- as.data.frame(
    matrix(
      rcpp_cx_report(bam.processed, pass, ctx), ncol=6, byrow=TRUE,
      dimnames=list(NULL, c("rname","strand","pos","context","meth","unmeth"))
    )
  )
  data.table::setDT(cx.report)
  cx.report[, data.table::setattr(rname,  "class", "factor")]
  cx.report[, data.table::setattr(rname,  "levels", levels(bam.processed$rname))]
  cx.report[, data.table::setattr(strand, "class", "factor")]
  cx.report[, data.table::setattr(strand, "levels", levels(bam.processed$strand))]
  cx.report[, `:=` (context = rcpp_char_to_context(context))]

  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(cx.report)
}


################################################################################

# descr: BED-assisted (amplicon/capture) report

.getBedReport <- function (bam.processed, pass, bed, bed.type,
                           match.tolerance, match.min.overlap,
                           verbose)
{
  if (verbose) message("Preparing ", bed.type, " report", appendLF=FALSE)
  tm <- proc.time()
  
  bam.subset <- data.table::data.table(
    strand=bam.processed$strand,
    bedmatch=.matchTarget(bam.processed=bam.processed, bed=bed,
                          bed.type=bed.type, match.tolerance=match.tolerance,
                          match.min.overlap=match.min.overlap),
    pass=factor(pass, levels=c(TRUE,FALSE))
  )
  data.table::setkey(bam.subset, bedmatch)
  bam.dt <- data.table::dcast(
    bam.subset[, list(nreads=.N), by=list(bedmatch,strand,pass), nomatch=0],
    bedmatch~pass+strand, value.var=c("nreads"), sep="", drop=FALSE, fill=0
  )
  bam.dt[,`:=` (`nreads+`=`FALSE+`+`TRUE+`,
                `nreads-`=`FALSE-`+`TRUE-`,
                VEF=(`TRUE+`+`TRUE-`)/(`FALSE+`+`TRUE+`+`FALSE-`+`TRUE-`) )]
  bed.dt <- data.table::as.data.table(bed)
  # bed.cl <- colnames(bed.dt)
  bed.dt[, bedmatch:=.I]
  bed.report <- data.table::merge.data.table(bed.dt, bam.dt, by="bedmatch",
                                             all=TRUE)[order(bedmatch)]
  
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(bed.report[,setdiff(names(bed.report),
                             c("bedmatch","FALSE+","FALSE-","TRUE+","TRUE-")),
                    with=FALSE])
}

################################################################################

# descr: calculates beta values and returns ECDF functions for BED file entries
# value: list of lists with context and out-of-context ECDF functions

.getBedEcdf <- function (bam.processed, bed, bed.type, bed.rows,
                         match.tolerance, match.min.overlap,
                         ctx.meth, ctx.unmeth, ooctx.meth, ooctx.unmeth,
                         verbose)
{
  if (verbose) message("Computing ECDFs for within- and out-of-context",
                       " per-read beta values", appendLF=FALSE)
  tm <- proc.time()
  
  bed.match <- .matchTarget(bam.processed=bam.processed, bed=bed,
                            bed.type=bed.type, match.tolerance=match.tolerance,
                            match.min.overlap=match.min.overlap)
  
  # Rcpp::sourceCpp("rcpp_get_xm_beta.cpp")
  ctx.beta=rcpp_get_xm_beta(bam.processed, ctx.meth, ctx.unmeth)
  ooctx.beta=rcpp_get_xm_beta(bam.processed, ooctx.meth, ooctx.unmeth)

  all.bed.rows <- sort(unique(bed.match), na.last=TRUE)
  if (is.null(bed.rows))
    bed.rows <- all.bed.rows
  else
    bed.rows <- intersect(bed.rows, all.bed.rows)
  
  bed.match.isna <- is.na(bed.match)
  bed.ecdf <- lapply(bed.rows, function (n) {
    matched.rows <- if (is.na(n))
      c(bed.match.isna)
    else
      c(!bed.match.isna & bed.match==n)
    return(c(context=stats::ecdf(ctx.beta[matched.rows]),
             out.of.context=stats::ecdf(ooctx.beta[matched.rows])))
  })
  names(bed.ecdf) <- as.character(as.character(bed)[bed.rows])
 
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(bed.ecdf)
}

################################################################################

# descr: calculates base frequences at particular positions
# value: data.table with base freqs

.getBaseFreqReport <- function (bam.processed, pass, vcf,
                                verbose)
{
  if (verbose) message("Extracting base frequences", appendLF=FALSE)
  tm <- proc.time()
  
  vcf.ranges <- BiocGenerics::sort(SummarizedExperiment::rowRanges(vcf))
  vcf.ranges <- vcf.ranges[BiocGenerics::width(vcf.ranges)==1 &
                           vapply(as.character(vcf.ranges$ALT),
                                  stringi::stri_length,
                                  FUN.VALUE=numeric(1), USE.NAMES=FALSE)==1]
  # GenomeInfoDb::seqlevels(vcf.ranges, pruning.mode="coarse") <-
  #   levels(bam.processed$rname)
  vcf.dt <- data.table::as.data.table(vcf.ranges)
  vcf.dt[, seqnames := factor(seqnames, levels=levels(bam.processed$rname))]
  
  freqs <- rcpp_get_base_freqs(bam.processed, pass, vcf.dt)
  colnames(freqs) <- c("","U+A","","U+C","U+T","","U+N","U+G",
                       "","U-A","","U-C","U-T","","U-N","U-G",
                       "","M+A","","M+C","M+T","","M+N","M+G",
                       "","M-A","","M-C","M-T","","M-N","M-G")
  
  bf.report <- data.table::data.table(
    name=names(vcf.ranges),
    vcf.dt[,.(seqnames, range=start, REF, ALT)],
    freqs[,grep("[ACTG]$",colnames(freqs))]
  )
  
  bf.report[REF=="A" & ALT=="C", `:=` (`M+Ref`=`M+A`,       `U+Ref`=`U+A`,       `M-Ref`=`M-A`,       `U-Ref`=`U-A`,
                                       `M+Alt`=`M+C`+`M+T`, `U+Alt`=`U+C`+`U+T`, `M-Alt`=`M-C`,       `U-Alt`=`U-C`      )]
  bf.report[REF=="A" & ALT=="T", `:=` (`M+Ref`=`M+A`,       `U+Ref`=`U+A`,       `M-Ref`=`M-A`,       `U-Ref`=`U-A`,
                                       `M+Alt`=`M+T`,       `U+Alt`=`U+T`,       `M-Alt`=`M-T`,       `U-Alt`=`U-T`      )]
  bf.report[REF=="A" & ALT=="G", `:=` (`M+Ref`=`M+A`,       `U+Ref`=`U+A`,       `M-Ref`=NA,          `U-Ref`=NA,
                                       `M+Alt`=`M+G`,       `U+Alt`=`U+G`,       `M-Alt`=NA,          `U-Alt`=NA         )]
  bf.report[REF=="C" & ALT=="A", `:=` (`M+Ref`=`M+C`+`M+T`, `U+Ref`=`U+C`+`U+T`, `M-Ref`=`M-C`,       `U-Ref`=`U-C`,
                                       `M+Alt`=`M+A`,       `U+Alt`=`U+A`,       `M-Alt`=`M-A`,       `U-Alt`=`U-A`      )]
  bf.report[REF=="C" & ALT=="T", `:=` (`M+Ref`=NA,          `U+Ref`=NA,          `M-Ref`=`M-C`,       `U-Ref`=`U-C`,
                                       `M+Alt`=NA,          `U+Alt`=NA,          `M-Alt`=`M-T`,       `U-Alt`=`U-T`      )]
  bf.report[REF=="C" & ALT=="G", `:=` (`M+Ref`=`M+C`+`M+T`, `U+Ref`=`U+C`+`U+T`, `M-Ref`=`M-C`,       `U-Ref`=`U-C`,
                                       `M+Alt`=`M+G`,       `U+Alt`=`U+G`,       `M-Alt`=`M-A`+`M-G`, `U-Alt`=`U-A`+`U-G`)]
  bf.report[REF=="T" & ALT=="A", `:=` (`M+Ref`=`M+T`,       `U+Ref`=`U+T`,       `M-Ref`=`M-T`,       `U-Ref`=`U-T`,
                                       `M+Alt`=`M+A`,       `U+Alt`=`U+A`,       `M-Alt`=`M-A`,       `U-Alt`=`U-A`      )]
  bf.report[REF=="T" & ALT=="C", `:=` (`M+Ref`=NA,          `U+Ref`=NA,          `M-Ref`=`M-T`,       `U-Ref`=`U-T`,
                                       `M+Alt`=NA,          `U+Alt`=NA,          `M-Alt`=`M-C`,       `U-Alt`=`U-C`      )]
  bf.report[REF=="T" & ALT=="G", `:=` (`M+Ref`=`M+T`,       `U+Ref`=`U+T`,       `M-Ref`=`M-T`,       `U-Ref`=`U-T`,
                                       `M+Alt`=`M+G`,       `U+Alt`=`U+G`,       `M-Alt`=`M-A`+`M-G`, `U-Alt`=`U-A`+`U-G`)]
  bf.report[REF=="G" & ALT=="A", `:=` (`M+Ref`=`M+G`,       `U+Ref`=`U+G`,       `M-Ref`=NA,          `U-Ref`=NA,
                                       `M+Alt`=`M+A`,       `U+Alt`=`U+A`,       `M-Alt`=NA,          `U-Alt`=NA         )]
  bf.report[REF=="G" & ALT=="C", `:=` (`M+Ref`=`M+G`,       `U+Ref`=`U+G`,       `M-Ref`=`M-A`+`M-G`, `U-Ref`=`U-A`+`U-G`,
                                       `M+Alt`=`M+C`+`M+T`, `U+Alt`=`U+C`+`U+T`, `M-Alt`=`M-C`,       `U-Alt`=`U-C`      )]
  bf.report[REF=="G" & ALT=="T", `:=` (`M+Ref`=`M+G`,       `U+Ref`=`U+G`,       `M-Ref`=`M-A`+`M-G`, `U-Ref`=`U-A`+`U-G`,
                                       `M+Alt`=`M+T`,       `U+Alt`=`U+T`,       `M-Alt`=`M-T`,       `U-Alt`=`U-T`      )]
  bf.report[, `:=` (SumRef=rowSums(bf.report[,.(`M+Ref`,`U+Ref`,`M-Ref`,`U-Ref`)], na.rm=TRUE),
                    SumAlt=rowSums(bf.report[,.(`M+Alt`,`U+Alt`,`M-Alt`,`U-Alt`)], na.rm=TRUE))]
  # FEp <- function (x) { if (any(is.na(x))) NA else stats::fisher.test(matrix(x, nrow=2))$p.value }
  # bf.report[, `:=` (`FEp+`=apply(bf.report[,.(`M+Ref`,`U+Ref`,`M+Alt`,`U+Alt`)], 1, FEp),
  #                   `FEp-`=apply(bf.report[,.(`M-Ref`,`U-Ref`,`M-Alt`,`U-Alt`)], 1, FEp))]
  bf.report[, `:=` (`FEp+`=rcpp_fep(bf.report, c("M+Ref","U+Ref","M+Alt","U+Alt")),
                    `FEp-`=rcpp_fep(bf.report, c("M-Ref","U-Ref","M-Alt","U-Alt")))]
  
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(bf.report)
}

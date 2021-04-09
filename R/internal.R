#' @importFrom Rsamtools testPairedEndBam 
#' @importFrom Rsamtools scanBamFlag
#' @importFrom Rsamtools bamFlagAND
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBam
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table as.data.table
#' @importFrom data.table dcast
#' @importFrom data.table merge.data.table
#' @importFrom data.table setorder
#' @importFrom data.table setkey
#' @importFrom stringi stri_length
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges reduce
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom VariantAnnotation ScanVcfParam
#' @importFrom VariantAnnotation readVcf
#' @importFrom VariantAnnotation expand
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @importFrom stats ecdf
#' @importFrom stats fisher.test
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

.onUnload <- function (libpath) {
  library.dynam.unload("epialleleR", libpath)
}

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

# descr: reads BAM files using Rsamtools
# value: unprocessed list output from Rsamtools::scanBam

.readBam <- function (bam.file,
                      min.mapq,
                      skip.duplicates,
                      verbose)
{
  if (verbose) message("Reading BAM file", appendLF=FALSE)
  tm <- proc.time()
    
  paired <- Rsamtools::testPairedEndBam(file=bam.file, index=NULL)
  
  flag <- Rsamtools::scanBamFlag()
  if (skip.duplicates)
    flag <- Rsamtools::bamFlagAND(flag, Rsamtools::scanBamFlag(isDuplicate=FALSE))
  if (paired)
    flag <- Rsamtools::bamFlagAND(flag, Rsamtools::scanBamFlag(isPaired=TRUE, isProperPair=TRUE))

  param <- Rsamtools::ScanBamParam(
    flag=flag,
    what=c("qname","flag","rname","strand","pos", "cigar", "seq",
           if (paired) c("mpos", "isize")),
    tag=c("XM","XG"),
    mapqFilter=min.mapq
  )
  
  bam <- Rsamtools::scanBam(bam.file, param=param)
  
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(bam)
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
  
  bed.df <- data.table::fread(file=bed.file, sep="\t", blank.lines.skip=TRUE, data.table=FALSE)
  colnames(bed.df)[1:3] <- c("chr", "start", "end")
  bed    <- GenomicRanges::makeGRangesFromDataFrame(bed.df, keep.extra.columns=TRUE, ignore.strand=TRUE,
                                                    starts.in.df.are.0based=zero.based.bed)
  
  
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
    vcf.genome <- stats::setNames(rep("unknown", length(bed.levels)), bed.levels)
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
  
  # if (gzip)
  #   report.file <- base::gzfile(base::sub("(\\.gz)?$", ".gz", report.file, ignore.case=TRUE), "w")
  # write.table(report, file=report.file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  data.table::fwrite(report, file=report.file, quote=FALSE, sep="\t",
                     row.names=FALSE, col.names=TRUE,
                     compress=if (gzip) "gzip" else "none")
  
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
}

################################################################################
# Functions: processing
################################################################################

# descr: process BAM data, merge reads if necessary
# value: data.table with fields qname, rname, strand, start, XM

.processBam <- function (bam,
                         verbose)
{
  if (verbose) message("Preprocessing BAM data:")
  tm <- proc.time()
  
  if (any(
    sapply(c(bam[[1]][c("qname","flag","rname","strand","pos","cigar","seq")],
             bam[[1]][[c("tag")]]), is.null)
  )) stop("BAM list object must contain data for the following BAM fields: ",
          "'qname', 'flag', 'rname', 'strand', 'pos', 'cigar', 'seq' ",
          "and the following tags: 'XM', 'XG'. ",
          "For paired-end sequencing files following BAM fields should be ",
          "present as well: 'mpos/pnext', 'isize/tlen'.")
  
  bam.data <- data.table::as.data.table(data.frame(bam, stringsAsFactors=FALSE))
  colnames(bam.data) <- sub("^tag.X", "X", colnames(bam.data))
  
  if (verbose) message("  Transforming sequences", appendLF=FALSE)
  bam.data[, `:=` (isfirst = bitwAnd(flag,128)==0,
                   seq     = rcpp_apply_cigar(cigar, seq),
                   XM      = rcpp_apply_cigar(cigar, XM))]
 
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)

  # fast merge reads, vectorised
  if (any(colnames(bam.data)=="mpos")) {
    if (verbose) message("  Merging pairs", appendLF=FALSE)
    tm <- proc.time()
    
    bam.data.first  <- bam.data[isfirst==TRUE, .(qname, rname, XG, pos, mpos, seq, XM, isize)]
    bam.data.second <- bam.data[isfirst!=TRUE, .(qname, seq, XM)]
    if (!identical(bam.data.first$qname, bam.data.second$qname))
      stop("Ungrouped reads? Please sort input BAM file by QNAME using 'samtools -n -o out.bam in.bam'")
    bam.data.first[, `:=` (rname  = factor(rname),
                           strand = factor(XG, levels=c("CT","GA"), labels=c("+","-")),
                           start  = base::pmin.int(pos, mpos),
                           width  = abs(isize),
                           seq    = rcpp_merge_ends(pos, seq, mpos, bam.data.second$seq, isize, 'N'),
                           XM     = rcpp_merge_ends(pos, XM,  mpos, bam.data.second$XM,  isize, '-') )]
    data.table::setorder(bam.data.first, rname, start)
    if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
    return(bam.data.first[,.(qname, rname, strand, start, seq, XM, width)])
  } else {
    bam.data[, width:=stringi::stri_length(XM)]
    data.table::setorder(bam.data, rname, pos)
    return(bam.data[,.(qname, rname, strand=XG, start=pos, seq, XM, width)])
  }
}

################################################################################

# descr: apply thresholding criteria to processed BAM reads
# value: bool vector with true for reads passing the threshold

.thresholdReads <- function (bam.processed,
                             ctx.meth, ctx.unmeth, ooctx.meth, ooctx.unmeth,
                             min.context.sites, min.context.beta, max.outofcontext.beta,
                             verbose)
{
  if (verbose) message("Thresholding reads", appendLF=FALSE)
  tm <- proc.time()
  
  # fast thresholding, vectorised
  pass <- rcpp_threshold_reads(
    bam.processed$XM,
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
  if (bed.type=="amplicon") {
    bed.match <- rcpp_match_amplicon(
      as.character(bam.processed$rname), bam.processed$start, bam.processed$start+bam.processed$width-1,
      as.character(GenomicRanges::seqnames(bed)), BiocGenerics::start(bed), BiocGenerics::end(bed),
      match.tolerance)
  } else if (bed.type=="capture") {
    bed.match <- rcpp_match_capture(
      as.character(bam.processed$rname), bam.processed$start, bam.processed$start+bam.processed$width-1,
      as.character(GenomicRanges::seqnames(bed)), BiocGenerics::start(bed), BiocGenerics::end(bed),
      match.min.overlap)
  }
  
  return(bed.match)
}

################################################################################
# Functions: reporting
################################################################################

# descr: prepare cytosine report for processed reads according to filter
# value: data.table with Bismark-formatted cytosine report

.getCytosineReport <- function (bam.processed,
                                ctx,
                                verbose)
{
  if (verbose) message("Preparing cytosine report", appendLF=FALSE)
  tm <- proc.time()
  
  # check if ordered? reorder if not
  # reporting, vectorised
  cx.report <- data.table::as.data.table(
    matrix(
      rcpp_cx_report(as.integer(bam.processed$rname), as.integer(bam.processed$strand),
                     bam.processed$start, bam.processed$XM, bam.processed$pass, ctx),
      ncol=6, byrow=TRUE, dimnames=list(NULL, c("rname","pos","strand","context","meth","unmeth"))
    )
  )
  cx.report[, `:=` (rname   = levels(bam.processed$rname )[rname ],
                    strand  = levels(bam.processed$strand)[strand],
                    context = rcpp_char_to_context(context),
                    triad   = "NNN")]

  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(cx.report[,.(rname, strand, pos, context, meth, unmeth, triad)])
}


################################################################################

# descr: BED-assisted (amplicon/capture) report

.getBedReport <- function (bam.processed, bed, bed.type,
                           match.tolerance, match.min.overlap,
                           verbose)
{
  if (verbose) message("Preparing ", bed.type, " report", appendLF=FALSE)
  tm <- proc.time()
  
  bam.processed[, `:=` (bedmatch=.matchTarget(bam.processed=bam.processed,
                                              bed=bed, bed.type=bed.type,
                                              match.tolerance=match.tolerance,
                                              match.min.overlap=match.min.overlap),
                        pass=factor(pass, levels=c(TRUE,FALSE)))]
  data.table::setkey(bam.processed, bedmatch)
  bam.dt <- data.table::dcast(bam.processed[, list(nreads=.N), by=list(bedmatch, strand, pass), nomatch=0],
                              bedmatch~pass+strand, value.var=c("nreads"), sep="", drop=FALSE, fill=0)
  bam.dt[,`:=` (`nreads+`=`FALSE+`+`TRUE+`,
                `nreads-`=`FALSE-`+`TRUE-`,
                VEF=(`TRUE+`+`TRUE-`)/(`FALSE+`+`TRUE+`+`FALSE-`+`TRUE-`) )]
  bed.dt <- data.table::as.data.table(bed)
  bed.cl <- colnames(bed.dt)
  bed.dt[, bedmatch:=.I]
  bed.report <- data.table::merge.data.table(bed.dt, bam.dt, by="bedmatch", all=TRUE)[order(bedmatch)]
  
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
  if (verbose) message("Computing ECDFs for within- and out-of-context per-read beta values", appendLF=FALSE)
  tm <- proc.time()
  
  bed.match <- .matchTarget(bam.processed=bam.processed, bed=bed, bed.type=bed.type,
                            match.tolerance=match.tolerance, match.min.overlap=match.min.overlap)
  
  # Rcpp::sourceCpp("rcpp_get_xm_beta.cpp")
  ctx.beta=rcpp_get_xm_beta(bam.processed$XM, ctx.meth, ctx.unmeth)
  ooctx.beta=rcpp_get_xm_beta(bam.processed$XM, ooctx.meth, ooctx.unmeth)

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

.getBaseFreqReport <- function (bam.processed, vcf,
                                verbose)
{
  if (verbose) message("Extracting base frequences", appendLF=FALSE)
  tm <- proc.time()
  
  vcf.ranges <- BiocGenerics::sort(SummarizedExperiment::rowRanges(vcf))
  vcf.ranges <- vcf.ranges[BiocGenerics::width(vcf.ranges)==1 &
                           sapply(as.character(vcf.ranges$ALT),
                                  stringi::stri_length, USE.NAMES=FALSE)==1]
  GenomeInfoDb::seqlevels(vcf.ranges, pruning.mode="coarse") <-
    levels(bam.processed$rname)
  freqs <- rcpp_get_base_freqs(as.integer(bam.processed$rname),
                               as.integer(bam.processed$strand),
                               bam.processed$start,
                               bam.processed$start+bam.processed$width-1,
                               bam.processed$seq,
                               bam.processed$pass,
                               as.integer(GenomicRanges::seqnames(vcf.ranges)),
                               BiocGenerics::start(vcf.ranges))
  colnames(freqs) <- c("","U+A","","U+C","U+T","","U+N","U+G",
                       "","U-A","","U-C","U-T","","U-N","U-G",
                       "","M+A","","M+C","M+T","","M+N","M+G",
                       "","M-A","","M-C","M-T","","M-N","M-G")
  
  bf.report <- data.table::data.table(
    name=names(vcf.ranges),
    seqnames=as.character(GenomicRanges::seqnames(vcf.ranges)),
    range=BiocGenerics::start(vcf.ranges),
    REF=as.character(vcf.ranges$REF),
    ALT=as.character(vcf.ranges$ALT),
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
  FEp <- function (x) { if (any(is.na(x))) NA else stats::fisher.test(matrix(x, nrow=2))$p.value }
  bf.report[, `:=` (`FEp+`=apply(bf.report[,.(`M+Ref`,`U+Ref`,`M+Alt`,`U+Alt`)], 1, FEp),
                    `FEp-`=apply(bf.report[,.(`M-Ref`,`U-Ref`,`M-Alt`,`U-Alt`)], 1, FEp))]
  
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(bf.report)
}


################################################################################
################################################################################
###
###
### Temporary copies below
###
###
################################################################################
################################################################################
# 
# 
# .getBaseFreqReport <- function (bam.processed, vcf,
#                                 verbose)
# {
#   if (verbose) message("Extracting base frequences", appendLF=FALSE)
#   tm <- proc.time()
#   
#   vcf.ranges <- unique(BiocGenerics::sort(SummarizedExperiment::rowRanges(vcf)))
#   vcf.ranges <- vcf.ranges[BiocGenerics::width(vcf.ranges)==1 &
#                              sapply(IRanges::CharacterList(vcf.ranges$ALT), length)==1]
#   vcf.ranges <- vcf.ranges[sapply(IRanges::CharacterList(vcf.ranges$ALT), stringi::stri_length)==1]
#   GenomeInfoDb::seqlevels(vcf.ranges) <- levels(bam.processed$rname)
#   freqs <- rcpp_get_base_freqs(as.integer(bam.processed$rname),
#                                as.integer(bam.processed$strand),
#                                bam.processed$start,
#                                bam.processed$start+bam.processed$width-1,
#                                bam.processed$seq,
#                                bam.processed$pass,
#                                as.integer(GenomicRanges::seqnames(vcf.ranges)),
#                                BiocGenerics::start(vcf.ranges))
#   colnames(freqs) <- c("","U+A","","U+C","U+T","","U+N","U+G",
#                        "","U-A","","U-C","U-T","","U-N","U-G",
#                        "","M+A","","M+C","M+T","","M+N","M+G",
#                        "","M-A","","M-C","M-T","","M-N","M-G")
#   
#   bf.report <- data.table::data.table(
#     name=names(vcf.ranges),
#     seqnames=as.character(GenomicRanges::seqnames(vcf.ranges)),
#     range=BiocGenerics::start(vcf.ranges),
#     REF=as.character(vcf.ranges$REF),
#     ALT=sapply(IRanges::CharacterList(vcf.ranges$ALT), paste, collapse=","),
#     freqs[,grep("[ACTG]$",colnames(freqs))]
#   )
#   
#   bf.report[REF=="A" & ALT=="C", `:=` (`M+Ref`=`M+A`,       `U+Ref`=`U+A`,       `M-Ref`=`M-A`,       `U-Ref`=`U-A`,
#                                        `M+Alt`=`M+C`+`M+T`, `U+Alt`=`U+C`+`U+T`, `M-Alt`=`M-C`,       `U-Alt`=`U-C`      )]
#   bf.report[REF=="A" & ALT=="T", `:=` (`M+Ref`=`M+A`,       `U+Ref`=`U+A`,       `M-Ref`=`M-A`,       `U-Ref`=`U-A`,
#                                        `M+Alt`=`M+T`,       `U+Alt`=`U+T`,       `M-Alt`=`M-T`,       `U-Alt`=`U-T`      )]
#   bf.report[REF=="A" & ALT=="G", `:=` (`M+Ref`=`M+A`,       `U+Ref`=`U+A`,       `M-Ref`=NA,          `U-Ref`=NA,
#                                        `M+Alt`=`M+G`,       `U+Alt`=`U+G`,       `M-Alt`=NA,          `U-Alt`=NA         )]
#   bf.report[REF=="C" & ALT=="A", `:=` (`M+Ref`=`M+C`+`M+T`, `U+Ref`=`U+C`+`U+T`, `M-Ref`=`M-C`,       `U-Ref`=`U-C`,
#                                        `M+Alt`=`M+A`,       `U+Alt`=`U+A`,       `M-Alt`=`M-A`,       `U-Alt`=`U-A`      )]
#   bf.report[REF=="C" & ALT=="T", `:=` (`M+Ref`=NA,          `U+Ref`=NA,          `M-Ref`=`M-C`,       `U-Ref`=`U-C`,
#                                        `M+Alt`=NA,          `U+Alt`=NA,          `M-Alt`=`M-T`,       `U-Alt`=`U-T`      )]
#   bf.report[REF=="C" & ALT=="G", `:=` (`M+Ref`=`M+C`+`M+T`, `U+Ref`=`U+C`+`U+T`, `M-Ref`=`M-C`,       `U-Ref`=`U-C`,
#                                        `M+Alt`=`M+G`,       `U+Alt`=`U+G`,       `M-Alt`=`M-A`+`M-G`, `U-Alt`=`U-A`+`U-G`)]
#   bf.report[REF=="T" & ALT=="A", `:=` (`M+Ref`=`M+T`,       `U+Ref`=`U+T`,       `M-Ref`=`M-T`,       `U-Ref`=`U-T`,
#                                        `M+Alt`=`M+A`,       `U+Alt`=`U+A`,       `M-Alt`=`M-A`,       `U-Alt`=`U-A`      )]
#   bf.report[REF=="T" & ALT=="C", `:=` (`M+Ref`=NA,          `U+Ref`=NA,          `M-Ref`=`M-T`,       `U-Ref`=`U-T`,
#                                        `M+Alt`=NA,          `U+Alt`=NA,          `M-Alt`=`M-C`,       `U-Alt`=`U-C`      )]
#   bf.report[REF=="T" & ALT=="G", `:=` (`M+Ref`=`M+T`,       `U+Ref`=`U+T`,       `M-Ref`=`M-T`,       `U-Ref`=`U-T`,
#                                        `M+Alt`=`M+G`,       `U+Alt`=`U+G`,       `M-Alt`=`M-A`+`M-G`, `U-Alt`=`U-A`+`U-G`)]
#   bf.report[REF=="G" & ALT=="A", `:=` (`M+Ref`=`M+G`,       `U+Ref`=`U+G`,       `M-Ref`=NA,          `U-Ref`=NA,
#                                        `M+Alt`=`M+A`,       `U+Alt`=`U+A`,       `M-Alt`=NA,          `U-Alt`=NA         )]
#   bf.report[REF=="G" & ALT=="C", `:=` (`M+Ref`=`M+G`,       `U+Ref`=`U+G`,       `M-Ref`=`M-A`+`M-G`, `U-Ref`=`U-A`+`U-G`,
#                                        `M+Alt`=`M+C`+`M+T`, `U+Alt`=`U+C`+`U+T`, `M-Alt`=`M-C`,       `U-Alt`=`U-C`      )]
#   bf.report[REF=="G" & ALT=="T", `:=` (`M+Ref`=`M+G`,       `U+Ref`=`U+G`,       `M-Ref`=`M-A`+`M-G`, `U-Ref`=`U-A`+`U-G`,
#                                        `M+Alt`=`M+T`,       `U+Alt`=`U+T`,       `M-Alt`=`M-T`,       `U-Alt`=`U-T`      )]
#   bf.report[, `:=` (SumRef=rowSums(bf.report[,.(`M+Ref`,`U+Ref`,`M-Ref`,`U-Ref`)], na.rm=TRUE),
#                     SumAlt=rowSums(bf.report[,.(`M+Alt`,`U+Alt`,`M-Alt`,`U-Alt`)], na.rm=TRUE))]
#   FEp <- function (x) { if (any(is.na(x))) NA else stats::fisher.test(matrix(x, nrow=2))$p.value }
#   bf.report[, `:=` (`FEp+`=apply(bf.report[,.(`M+Ref`,`U+Ref`,`M+Alt`,`U+Alt`)], 1, FEp),
#                     `FEp-`=apply(bf.report[,.(`M-Ref`,`U-Ref`,`M-Alt`,`U-Alt`)], 1, FEp))]
#   
#   if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
#   return(bf.report)
# }

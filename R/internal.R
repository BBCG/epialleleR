#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table data.table
#' @importFrom data.table as.data.table
#' @importFrom data.table dcast
#' @importFrom data.table merge.data.table
#' @importFrom data.table setorder
#' @importFrom data.table setkey
#' @importFrom data.table setDT
#' @importFrom data.table setattr
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
#' @importFrom utils packageVersion
#' @importFrom Rcpp evalCpp
#' @useDynLib epialleleR, .registration=TRUE


# internal globals, constants and helper functions 
#

################################################################################
# Globals, unload
################################################################################

utils::globalVariables(
  c(".", ".I", ".N", ":=", "bedmatch", "context", "rname", "start", "strand",
    "templid", "FALSE+", "FALSE-", "TRUE+", "TRUE-", "REF", "ALT",
    "M+Ref","U+Ref","M+Alt","U+Alt", "M-Ref","U-Ref","M-Alt","U-Alt",
    "M+A", "M+C", "M+G", "M+T", "M-A", "M-C", "M-G", "M-T",
    "U+A", "U+C", "U+G", "U+T", "U-A", "U-C", "U-G", "U-T")
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

# descr: check BAM file before reading (using HTSlib)
# value: TRUE (for paired-end) or FALSE (for single-end) or error
#        if BAM loading is not possible

.checkBam <- function (bam.file,
                       verbose)
{
  if (verbose) message("Checking BAM file: ", appendLF=FALSE)
  
  bam.file <- path.expand(bam.file)
  check.out <- rcpp_check_bam(bam.file)
  
  # main logic
  if (check.out$nrecs==0) {                                     # no records
    stop("Empty file provided! Exiting",
         call.=FALSE)
  } else if (is.null(check.out$XG) & !is.null(check.out$YD)) {  # YDs but no XGs
    stop("No XG tags found (though YD tags are there)! BWA-meth alignment?\n",
         "If so, make methylation calls using epialleleR::callMethylation.\n",
         "Exiting", call.=FALSE)
  } else if (is.null(check.out$XG) & !is.null(check.out$ZS)) {  # ZSs but no XGs
    stop("No XG tags found (though ZS tags are there)! BSMAP alignment?\n",
         "If so, make methylation calls using epialleleR::callMethylation.\n",
         "Exiting", call.=FALSE)
  } else if (is.null(check.out$XM) & !is.null(check.out$XG)) {  # XGs but no XMs
    stop("No XM tags found! Was methylation called successfully?\n",
         "If not, make methylation calls using epialleleR::callMethylation.\n",
         "Exiting", call.=FALSE)
  } else if (check.out$npaired < check.out$nrecs/2) {       # predominantly SE
    paired <- FALSE
  } else {
    if (check.out$ntempls*2 < check.out$npaired - 1) {      # not sorted by name
      stop("BAM file seems to be paired-end but not sorted by name!\n",
           "Please sort using 'samtools sort -n -o out.bam in.bam'.\n",
           "Exiting", call.=FALSE)
    }
    paired <- TRUE
  }
  
  if (verbose) message(ifelse(paired, "paired-end, name-sorted", "single-end"),
                       " alignment detected", appendLF=TRUE)
  return(paired)
}

################################################################################

# descr: reads genomic (bgzipped) FASTA files using HTSlib
# value: list

.readGenome <- function (genome.file,
                         nthreads,
                         verbose)
{
  if (verbose) message("Reading reference genome file ", appendLF=FALSE)
  tm <- proc.time()
  
  genome.file <- path.expand(genome.file)
  genome.processed <- rcpp_read_genome(genome.file, nthreads)
  
  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(genome.processed)
}

################################################################################

# descr: reads BAM files using HTSlib
# value: data.table

.readBam <- function (bam.file,
                      paired,
                      min.mapq,
                      min.baseq,
                      skip.duplicates,
                      nthreads,
                      verbose)
{
  if (verbose) message("Reading ", ifelse(paired, "paired", "single"), 
                       "-end BAM file ", appendLF=FALSE)
  tm <- proc.time()
  
  bam.file <- path.expand(bam.file)
  if (paired) {
    bam.processed <- rcpp_read_bam_paired(bam.file, min.mapq, min.baseq, 
                                          skip.duplicates, nthreads)
  } else {
    bam.processed <- rcpp_read_bam_single(bam.file, min.mapq, min.baseq, 
                                          skip.duplicates, nthreads)
  }
  
  data.table::setDT(bam.processed)
  bam.processed[,templid:=c(0:(.N-1))]
  data.table::setorder(bam.processed, rname, start)
  
  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(bam.processed)
}

################################################################################

# descr: (fast) reads BED file with amplicons
# value: object of type GRanges

.readBed <- function (bed.file,
                      zero.based.bed,
                      verbose)
{
  if (verbose) message("Reading BED file ", appendLF=FALSE)
  tm <- proc.time()
  
  bed.df <- data.table::fread(file=bed.file, sep="\t", blank.lines.skip=TRUE,
                              data.table=FALSE)
  colnames(bed.df)[seq_len(3)] <- c("chr", "start", "end")
  bed <- GenomicRanges::makeGRangesFromDataFrame(
    bed.df, keep.extra.columns=TRUE, ignore.strand=TRUE,
    starts.in.df.are.0based=zero.based.bed
  )
  
  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
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
  if (verbose) message("Reading VCF file ", appendLF=FALSE)
  tm <- proc.time()
  
  vcf.genome <- "unknown"
  vcf.param  <- VariantAnnotation::ScanVcfParam(info=NA, geno=NA)
  
  if (!is.null(bed)) {
    bed <- GenomicRanges::reduce(bed)
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
  
  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
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
  if (verbose) message("Writing the report ", appendLF=FALSE)
  tm <- proc.time()
  
  data.table::fwrite(report, file=report.file, quote=FALSE, sep="\t",
                     row.names=FALSE, col.names=TRUE,
                     compress=if (gzip) "gzip" else "none")
  
  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
}

################################################################################
# Functions: processing
################################################################################

# descr: writing out example BAM
# value: number of records written

.simulateBam <- function (output.bam.file,
                          qname, flag, rname,
                          pos, mapq, cigar, rnext,
                          pnext, tlen, seq, qual, ...,
                          verbose)
{
  if (verbose) message("Writing sample BAM ", appendLF=FALSE)
  tm <- proc.time()
  
  if (!is.null(output.bam.file)) output.bam.file <- path.expand(output.bam.file)
  tags <- list(...)
  nrecs <- max(
    sapply(c(as.list(environment()), tags), length)[names(match.call())],
    1, na.rm=TRUE
  )
  
  if (is.null(qname)) qname <- sprintf("q%.04i", seq_len(nrecs))
  else                qname <- rep_len(qname, nrecs)
  if (is.null(flag))  flag  <- rep_len(0, nrecs)
  else                flag  <- rep_len(flag, nrecs)
  if (is.null(rname)) rname <- factor(rep_len("chrS", nrecs))
  else                rname <- factor(rep_len(rname, nrecs))
  if (is.null(pos))   pos   <- rep_len(1, nrecs)
  else                pos   <- rep_len(pos, nrecs)
  if (is.null(mapq))  mapq  <- rep_len(60, nrecs)
  else                mapq  <- rep_len(mapq, nrecs)
  if (is.null(seq)) {
    if ("XM" %in% names(tags)) {nbases <- nchar(tags$XM)}
    else if (!is.null(tlen))   {nbases <- tlen}
    else                       {nbases <- 10}
    seq <- sapply(nbases, function (l) {
      paste(sample(c("A","C","T","G"), l, replace=TRUE), collapse="")
    })
  }
                      seq   <- rep_len(seq, nrecs)
  if (is.null(cigar)) cigar <- paste0(nchar(seq), "M")
  else                cigar <- rep_len(cigar, nrecs)
  if (is.null(rnext)) rnext <- factor(rep_len("chrS", nrecs))
  else                rnext <- factor(rep_len(rnext, nrecs))
  if (is.null(pnext)) pnext <- rep_len(1, nrecs)
  else                pnext <- rep_len(pnext, nrecs)
  if (is.null(tlen))  tlen  <- nchar(seq)
  else                tlen  <- rep_len(tlen, nrecs)
  if (is.null(qual))  qual  <- sapply(lapply(nchar(seq), rep_len, x="F"), paste, collapse="")
  else                qual  <- rep_len(qual, nrecs)
  
  header <- c(
    sprintf("@SQ\tSN:%s\tLN:%i", levels(rname), max(pos, pnext)+max(tlen)-1),
    sprintf("@PG\tID:epialleleR\tPN:epialleleR\tVN:%s\tCL:rcpp_simulate_bam()",
            utils::packageVersion("epialleleR"))
  )
  fields <- data.table::data.table(
    qname=qname, flag=flag, tid=as.integer(rname)-1, pos=pos-1, mapq=mapq,
    cigar=cigar, mtid=as.integer(rnext)-1, mpos=pnext-1,
    isize=tlen, seq=seq, qual=qual
  )
  i_tags <- data.table::as.data.table(tags[sapply(tags, is.integer)])
  i_tags <- i_tags[rep_len(seq_len(nrow(i_tags)), nrecs)]
  s_tags <- data.table::as.data.table(tags[sapply(tags, is.character)])
  s_tags <- s_tags[rep_len(seq_len(nrow(s_tags)), nrecs)]
  
  if (!is.null(output.bam.file))
    result <- rcpp_simulate_bam(header, fields, i_tags, s_tags, output.bam.file)
  else
    result <- cbind(fields, i_tags, s_tags)
  
  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(result)
}

################################################################################

# descr: calling methylation (writing XG/XM)
# value: simple statistics and output BAM

.callMethylation <- function (input.bam.file, output.bam.file,
                              genome, nthreads, verbose)
{
  if (verbose) message("Making methylation calls ", appendLF=FALSE)
  tm <- proc.time()
  
  input.bam.file <- path.expand(input.bam.file)
  output.bam.file <- path.expand(output.bam.file)
  check.out <- rcpp_check_bam(input.bam.file)
  
  if (check.out$nrecs==0) {                             # no records
    stop("Empty file provided! Exiting", call.=FALSE)
  } else if (!is.null(check.out$XG)) {
    tag <- "XG"
  } else if (!is.null(check.out$YD)) {
    tag <- "YD"
  } else if (!is.null(check.out$ZS)) {
    tag <- "ZS"
  } else stop("Unable to call methylation: neither of XG/YD/ZS tags is present",
              " (genome strand unknown).\nExiting", call.=FALSE)
  
  result <- rcpp_call_methylation_genome(
    input.bam.file, output.bam.file, genome, tag, nthreads
  )
  
  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(result)
}

################################################################################

# descr: apply thresholding criteria to processed BAM reads
# value: bool vector with true for reads passing the threshold

.thresholdReads <- function (bam.processed,
                             ctx.meth, ctx.unmeth, ooctx.meth, ooctx.unmeth,
                             min.context.sites, min.context.beta,
                             max.outofcontext.beta, verbose)
{
  if (verbose) message("Thresholding reads ", appendLF=FALSE)
  tm <- proc.time()
  
  # fast thresholding, vectorised
  pass <- rcpp_threshold_reads(
    bam.processed,
    ctx.meth, ctx.unmeth, ooctx.meth, ooctx.unmeth,
    min.context.sites, min.context.beta, max.outofcontext.beta
  )
  
  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
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
  if (verbose) message("Preparing cytosine report ", appendLF=FALSE)
  tm <- proc.time()
  
  # must be ordered
  cx.report <- rcpp_cx_report(bam.processed, pass, ctx)
  data.table::setDT(cx.report)
  
  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(cx.report)
}


################################################################################

# descr: prepare lMHL report for processed reads
# value: data.table with lMHL report

.getMhlReport <- function (bam.processed,
                           ctx, max.window, min.length,
                           verbose)
{
  if (verbose) message("Preparing lMHL report ", appendLF=FALSE)
  tm <- proc.time()
  
  # must be ordered
  mhl.report <- rcpp_mhl_report(bam.processed, ctx, max.window, min.length)
  data.table::setDT(mhl.report)
  
  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(mhl.report)
}


################################################################################

# descr: BED-assisted (amplicon/capture) report

.getBedReport <- function (bam.processed, pass, bed, bed.type,
                           match.tolerance, match.min.overlap,
                           verbose)
{
  if (verbose) message("Preparing ", bed.type, " report ", appendLF=FALSE)
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
  
  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
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
                       " per-read beta values ", appendLF=FALSE)
  tm <- proc.time()
  
  bed.match <- .matchTarget(bam.processed=bam.processed, bed=bed,
                            bed.type=bed.type, match.tolerance=match.tolerance,
                            match.min.overlap=match.min.overlap)
  
  # Rcpp::sourceCpp("rcpp_get_xm_beta.cpp")
  ctx.beta <- rcpp_get_xm_beta(bam.processed, ctx.meth, ctx.unmeth)
  ooctx.beta <- rcpp_get_xm_beta(bam.processed, ooctx.meth, ooctx.unmeth)

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
 
  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(bed.ecdf)
}

################################################################################

# descr: calculates base frequences at particular positions
# value: data.table with base freqs

.getBaseFreqReport <- function (bam.processed, pass, vcf,
                                verbose)
{
  if (verbose) message("Extracting base frequences ", appendLF=FALSE)
  tm <- proc.time()
  
  vcf.ranges <- BiocGenerics::sort(SummarizedExperiment::rowRanges(vcf))
  vcf.ranges <- vcf.ranges[BiocGenerics::width(vcf.ranges)==1 &
                           vapply(as.character(vcf.ranges$ALT),
                                  nchar,
                                  FUN.VALUE=numeric(1), USE.NAMES=FALSE)==1]
  # GenomeInfoDb::seqlevels(vcf.ranges, pruning.mode="coarse") <-
  #   levels(bam.processed$rname)
  vcf.dt <- data.table::as.data.table(vcf.ranges)
  vcf.dt[, seqnames := factor(seqnames, levels=levels(bam.processed$rname))]
  if (all(is.na(vcf.dt$seqnames)))
    stop("Looks like seqlevels styles of BAM and VCF don't match. ",
         "Please provide VCF as an object with correct seqlevels.")
  
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
  
  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(bf.report)
}

################################################################################

# descr: extracts methylation patterns for particular range
# value: data.table with patterns

.getPatterns <- function (bam.processed, bed, bed.row, match.min.overlap,
                          extract.context, min.context.freq,
                          clip.patterns, strand.offset, highlight.positions,
                          verbose)
{
  if (verbose) message("Extracting methylation patterns ", appendLF=FALSE)
  tm <- proc.time()
  
  bed.dt <- data.table::as.data.table(bed)[bed.row]
  bed.dt[, seqnames := factor(seqnames, levels=levels(bam.processed$rname))]
  
  highlight.positions <- sort(unique(
    highlight.positions[highlight.positions>=bed.dt$start &
                        highlight.positions<=bed.dt$end]
  ))
  
  patterns <- rcpp_extract_patterns(bam.processed,
                                    as.integer(bed.dt$seqnames),
                                    as.integer(bed.dt$start),
                                    as.integer(bed.dt$end),
                                    match.min.overlap, extract.context,
                                    min.context.freq,
                                    clip.patterns, strand.offset,
                                    highlight.positions)
  data.table::setDT(patterns)
  colnames(patterns) <- sub("^X([0-9]+)$", "\\1", colnames(patterns))
  
  if (verbose) message(sprintf("[%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(patterns)
}

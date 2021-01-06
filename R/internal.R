# internal constants and helper functions 
#

################################################################################
# Constants
################################################################################

# descr: methylation context to XM bases
#' noRd
.context.to.bases <- list (
  "CG"  = list (ctx.meth   = "Z",    ctx.unmeth   = "z",
                ooctx.meth = "XHU",  ooctx.unmeth = "xhu"),
  "CHG" = list (ctx.meth   = "X",    ctx.unmeth   = "x",
                ooctx.meth = "ZHU",  ooctx.unmeth = "zhu"),
  "CHH" = list (ctx.meth   = "H",    ctx.unmeth   = "h",
                ooctx.meth = "ZXU",  ooctx.unmeth = "zxu"),
  "CxG" = list (ctx.meth   = "ZX",   ctx.unmeth   = "zx",
                ooctx.meth = "HU",   ooctx.unmeth = "hu"),
  "CX"  = list (ctx.meth   = "ZXHU", ctx.unmeth   = "zxhu",
                ooctx.meth = "",     ooctx.unmeth = "")
)

################################################################################
# Functions
################################################################################

# descr: reads BAM files using Rsamtools
# value: unprocessed list output from Rsamtools::scanBam
#' noRd
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
    what=c("qname","flag","rname","strand","pos", "mapq","cigar",
           if (paired) c("mpos", "isize")),
    tag=c("XM","XG"),
    mapqFilter=min.mapq
  )
  
  bam <- Rsamtools::scanBam(bam.file, param=param)
  
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(bam)
}

################################################################################

# descr: process BAM data, merge reads if necessary
# value: tibble with fields qname, rname, strand, start, XM
#' noRd
.processBam <- function (bam,
                         verbose)
{
  if (verbose) message("Processing BAM data:")
  tm <- proc.time()
  
  if (any(
    sapply(c(bam[[1]][c("qname","flag","rname","strand","pos", "mapq","cigar")],
             bam[[1]][[c("tag")]]), is.null)
  )) stop("BAM list object must contain data for the following BAM fields: ",
          "'qname', 'flag', 'rname', 'strand', 'pos', 'mapq', 'cigar' ",
          "and following tags: 'XM', 'XG'. ",
          "For paired-end sequencing files following BAM fields should be ",
          "included as well: 'mpos', 'isize'.")
  
  bam.data <- dplyr::as_tibble(data.frame(bam, stringsAsFactors=FALSE),
                               .name_repair="minimal")
  colnames(bam.data) <- gsub("tag.", "", colnames(bam.data), fixed=TRUE)
  
  if (verbose) message("  Applying CIGARs to methylation tags", appendLF=FALSE)
  bam.data <- bam.data %>% 
    dplyr::mutate(
      rname=as.factor(rname),
      XG=factor(XG, levels=c("CT","GA"), labels=c("+","-")),
      isfirst=bitwAnd(flag,128)==0,
      XM.norm=as.character(
        GenomicAlignments::sequenceLayer(Biostrings::BStringSet(XM, use.names=FALSE), cigar),
        use.names=FALSE
      )
    )
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)

  # fast merge reads, vectorised
  # Rcpp::sourceCpp("rcpp_merge_ends.cpp")
  system.time(
    if (any(colnames(bam.data)=="mpos")) {
      if (verbose) message("  Merging read pairs", appendLF=FALSE)
      tm <- proc.time()
      
      bam.data.first  <- bam.data[bam.data$isfirst,]
      bam.data.second <- bam.data[!bam.data$isfirst,]
      if (!identical(bam.data.first$qname, bam.data.second$qname))
        stop("Ungrouped reads?? Please sort input BAM file by QNAME using 'samtools -n'")
      bam.data.merged <- bam.data.first %>%
        dplyr::mutate(start=base::pmin.int(pos, mpos),
                      width=abs(isize),
                      XM=rcpp_merge_ends_vector(bam.data.first$pos, bam.data.first$XM.norm,
                                                bam.data.second$pos, bam.data.second$XM.norm,
                                                bam.data.second$isize)) %>%
        dplyr::select(qname, rname, strand=XG, start, XM, width)
      
      if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
    } else {
      bam.data.merged <- bam.data %>%
        dplyr::mutate(width=stringi::stri_length(XM.norm)) %>%
        dplyr::select(qname, rname, strand=XG, start=pos, XM=XM.norm, width)
    }
  )
  
  return(bam.data.merged)
}

################################################################################

# descr: apply thresholding criteria to processed BAM reads
# value: bool vector with true for reads passing the threshold
#' noRd
.thresholdReads <- function (bam.processed,
                             ctx.meth, ctx.unmeth, ooctx.meth, ooctx.unmeth,
                             min.context.sites, min.context.beta, max.outofcontext.beta,
                             verbose)
{
  if (verbose) message("Thresholding reads", appendLF=FALSE)
  tm <- proc.time()
  
  # fast thresholding, vectorised
  # Rcpp::sourceCpp("rcpp_threshold_reads.cpp")
  pass <- rcpp_threshold_reads_vector(
    bam.processed$XM,
    ctx.meth, ctx.unmeth, ooctx.meth, ooctx.unmeth,
    min.context.sites, min.context.beta, max.outofcontext.beta
  )
  
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(pass)
}

################################################################################

# descr: prepare cytosine report for processed reads according to filter
# value: tibble with Bismark-formatted cytosine report
# TODO: rewrite "summarise" in C++
#   std::map<int, std::map<int, std::array<int,2>>>
#' noRd
.getCytosineReport <- function (bam.processed,
                                ctx,
                                verbose)
{
  if (verbose) message("Preparing cytosine report", appendLF=FALSE)
  tm <- proc.time()
  
  # reporting, vectorised
  # Rcpp::sourceCpp("rcpp_parse_xm.cpp")
  cx.report <- dplyr::as_tibble(matrix(
    rcpp_parse_xm_vector(as.integer(bam.processed$rname),
                         as.integer(bam.processed$strand),
                         bam.processed$start, bam.processed$XM,
                         bam.processed$pass, ctx),
      ncol=6, byrow=TRUE, dimnames=list(NULL, c("rname","strand","pos","ctx","meth","unmeth"))
    )) %>%
      dplyr::group_by(rname,pos,strand) %>%
      dplyr::summarise(meth=sum(meth),
                       unmeth=sum(unmeth),
                       ctx=ctx[1],
                       .groups="drop") %>%
      dplyr::mutate(rname=levels(bam.processed$rname)[rname],
                    strand=levels(bam.processed$strand)[strand],
                    ctx=rcpp_char_to_context_vector(ctx),
                    triad="NNN")
  
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(cx.report)
}

################################################################################

# descr: (fast) writes the report
# value: void
#' noRd
.writeReport <- function (report,
                          report.file,
                          gzip,
                          verbose)
{
  if (verbose) message("Writing the report", appendLF=FALSE)
  tm <- proc.time()
  
  if (gzip)
    report.file <- base::gzfile(base::sub("(\\.gz)?$", ".gz", report.file, ignore.case=TRUE), "w")
  write.table(report, file=report.file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
}

################################################################################

# descr: (fast) reads BED file with amplicons
# value: GRanges
#' noRd
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

# descr: matching BED target (amplicon/capture)
# value: numeric vector
#' noRd
.matchTarget <- function (bam.processed, bed, bed.type,
                          match.tolerance, match.min.overlap)
{
  # fast, vectorised
  # Rcpp::sourceCpp("rcpp_match_target.cpp")
  if (bed.type=="amplicon") {
    bed.match <- rcpp_match_amplicon_vector(
      as.character(bam.processed$rname), bam.processed$start, bam.processed$start+bam.processed$width-1,
      as.character(GenomicRanges::seqnames(bed)), GenomicRanges::start(bed), GenomicRanges::end(bed),
      match.tolerance)
  } else if (bed.type=="capture") {
    bed.match <- rcpp_match_capture_vector(
      as.character(bam.processed$rname), bam.processed$start, bam.processed$start+bam.processed$width-1,
      as.character(GenomicRanges::seqnames(bed)), GenomicRanges::start(bed), GenomicRanges::end(bed),
      match.min.overlap)
  }

  return(bed.match)
}

################################################################################

# descr: BED-assisted (amplicon/capture) report
# value: tibble
#' noRd
.getBedReport <- function (bam.processed, bed, bed.type,
                           match.tolerance, match.min.overlap,
                           verbose)
{
  if (verbose) message("Preparing ", bed.type, " report", appendLF=FALSE)
  tm <- proc.time()

  bed.match <- .matchTarget(bam.processed=bam.processed, bed=bed, bed.type=bed.type,
                            match.tolerance=match.tolerance, match.min.overlap=match.min.overlap)
  
  bed.report <- bam.processed %>%
    dplyr::mutate(pass=factor(pass, levels=c(TRUE,FALSE), labels=c("above","below")),
                  match=bed.match) %>%
    dplyr::group_by(match, pass, strand, .drop=FALSE) %>%
    dplyr::summarise(nreads=n(), .groups="drop") %>%
    data.table::as.data.table() %>%
    data.table::dcast(match~pass+strand, value.var=c("nreads"), sep="", drop=FALSE)
  
  bed.report <- dplyr::bind_cols(data.frame(bed, stringsAsFactors=FALSE)[bed.report$match,],
                                 #bed.report[,-c("match")],
                                 `nreads+`=bed.report$`above+` + bed.report$`below+`,
                                 `nreads-`=bed.report$`above-` + bed.report$`below-`,
                                 VEF=(bed.report$`above+` + bed.report$`above-`)/(bed.report$`above+` + bed.report$`below+` +
                                                                                  bed.report$`above-` + bed.report$`below-`)) %>%
    dplyr::as_tibble(.name_repair="minimal")

  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(bed.report)
}

################################################################################

# descr: calculates beta values and returns ECDF functions for BED file entries
# value: list of lists with context and out-of-context ECDF functions
#' noRd
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
  ctx.beta=rcpp_get_xm_beta_vector(bam.processed$XM, ctx.meth, ctx.unmeth)
  ooctx.beta=rcpp_get_xm_beta_vector(bam.processed$XM, ooctx.meth, ooctx.unmeth)

  all.bed.rows <- sort(unique(bed.match), na.last=FALSE)
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
    return(c(context=ecdf(ctx.beta[matched.rows]),
             out.of.context=ecdf(ooctx.beta[matched.rows])))
  })
  names(bed.ecdf) <- as.character(as.character(bed)[bed.rows])
 
  if (verbose) message(sprintf(" [%.3fs]",(proc.time()-tm)[3]), appendLF=TRUE)
  return(bed.ecdf)
}

################################################################################
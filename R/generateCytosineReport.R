#' generateCytosineReport
#'
#' @description
#' `generateCytosineReport` desc. Epiallele-aware
#'
#' @details
#' details
#'
#' @param bam param desc
#' @return return desc
#' @seealso whatelse
#' @examples
#'   sessionInfo
# @import parallel
# @importFrom Rsamtools 
#' @export
generateCytosineReport <- function (bam,
                                    report.file=NULL,
                                    threshold.reads=TRUE,
                                    threshold.context=c("CG", "CHG", "CHH", "CxG", "CX"),
                                    min.context.sites=2,
                                    min.context.beta=0.5,
                                    max.outofcontext.beta=0.1, # double 0 to 1
                                    report.context=threshold.context,
                                    min.mapq=0,
                                    skip.duplicates=FALSE, # http://www.htslib.org/doc/samtools-markdup.html https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/DuplicateMarking_fDG.htm
                                    gzip=FALSE,
                                    verbose=TRUE)
{
  # del afterwards ----------------------------------------------
    # library(Rsamtools)
    # library(GenomicAlignments)
    # library(Biostrings)
    library(dplyr)
    # library(Rcpp)
    report.file <- NULL
    min.mapq <- 0 #33
    skip.duplicates <- FALSE
    threshold.reads <- TRUE
    threshold.context <- "CG"
    min.context.sites <- 2
    min.context.beta <- 0.5
    max.outofcontext.beta <- 0.01
    report.context <- threshold.context
    gzip <- FALSE
    verbose <- TRUE
    bam <- "/scratch/epialleleR/w/WHIP050-D701-D501_S1_L001/WHIP050-D701-D501_S1_L001.bam"
    # bam <- "/scratch/epialleleR/b/WHIP050-D701-D501_S1_L001/WHIP050-D701-D501_S1_L001_pe.bam"
    # bam <- "/scratch/epialleleR/b/1006-PT2-M_S3_L001/1006-PT2-M_S3_L001_pe.bam"
    report.file <- base::sub("(\\.bam)?$", paste0(".",report.context,"_report.txt"), bam, ignore.case=TRUE)
    
    Rcpp::sourceCpp("/Home/siv22/oni062/work/packages/epialleleR/src/rcpp_merge_ends.cpp")
    Rcpp::sourceCpp("/Home/siv22/oni062/work/packages/epialleleR/src/rcpp_threshold_reads.cpp")
    Rcpp::sourceCpp("/Home/siv22/oni062/work/packages/epialleleR/src/rcpp_parse_xm.cpp")
    Rcpp::sourceCpp("/Home/siv22/oni062/work/packages/epialleleR/src/rcpp_char_to_context.cpp")
  # del afterwards ----------------------------------------------
  
  threshold.context <- match.arg(threshold.context, threshold.context)
  report.context    <- match.arg(report.context, report.context)
  
  if (is.character(bam))
    bam <- .readBam(bam.file=bam, min.mapq=min.mapq,
                    skip.duplicates=skip.duplicates, verbose=verbose)
  
  bam.processed <- .processBam(bam=bam, verbose=verbose)
  
  if (threshold.reads) {
    bam.processed$pass <- .thresholdReads(
      bam.processed=bam.processed,
      ctx.meth=.context.to.bases[[threshold.context]][["ctx.meth"]],
      ctx.unmeth=.context.to.bases[[threshold.context]][["ctx.unmeth"]],
      ooctx.meth=.context.to.bases[[threshold.context]][["ooctx.meth"]],
      ooctx.unmeth=.context.to.bases[[threshold.context]][["ooctx.unmeth"]],
      min.context.sites=min.context.sites,
      min.context.beta=min.context.beta,
      max.outofcontext.beta=max.outofcontext.beta,
      verbose=verbose
    )
  } else {
    bam.processed$pass <- TRUE
  }
  
  cx.report <- .getCytosineReport(
    bam.processed=bam.processed,
    ctx=paste0(.context.to.bases[[report.context]][["ctx.meth"]],
               .context.to.bases[[report.context]][["ctx.unmeth"]]),
    verbose=verbose
  )
  
  if (is.null(report.file))
    return(cx.report)
  else
    .writeReport(report=cx.report, report.file=report.file, gzip=gzip,
                 verbose=verbose)
}


#|[c]{}^

# ##############################################################################

# if (!identical(bam.data[ bam.data$isfirst,"qname"],
#                bam.data[!bam.data$isfirst,"qname"])) {
#   warning("Unsorted bam. Sorting...")
#   bam.data <- bam.data[order(bam.data$qname, bam.data$isfirst),]
# }

# # normalizing methylation tags according to CIGAR
# system.time(
#   bam.data$XM.norm <- as.character(GenomicAlignments::sequenceLayer(Biostrings::BStringSet(bam.data$XM), bam.data$cigar))
# )

# # slow attemp to join them
# joinSeq <- function (isfirst, basepos, width, seq) {
#   joined  <- stringi::stri_dup(".",width[1])
#   substring(joined, basepos[!isfirst]) <- seq[!isfirst]
#   substring(joined, basepos[ isfirst]) <- seq[ isfirst]
#   return(joined)
#   # return(stringi::stri_flatten(seq))
# }
# system.time(
#   bam.full <- bam.data %>% 
#     dplyr::group_by(qname, rname, XG, start, width) %>%
#     dplyr::summarise(XM=joinSeq(isfirst, pos-start+1, width, XM.norm),
#                      .groups="drop")
# )
# bam.full.h <- head(bam.full)

# # truncating methylation tag of 2nd read
# system.time(
#   bam.data.first <- bam.data[bam.data$isfirst,] %>%
#     dplyr::mutate(start=pos)
# )
# rownames(bam.data.first) <- bam.data.first$qname
# system.time(
#   bam.data.second.leftmost  <- bam.data[!bam.data$isfirst & bam.data$isize>0,] %>%
#     dplyr::mutate(start=pos,
#                   XM.norm=stringi::stri_sub(XM.norm, to=mpos-pos))
# )
# # rownames(bam.data.second.leftmost) <- bam.data.second.leftmost$qname
# system.time(
#   bam.data.second.rightmost <- bam.data[!bam.data$isfirst & bam.data$isize<0,] %>%
#     dplyr::mutate(start=base::pmax.int(pos, mpos + nchar(bam.data.first[qname,"XM.norm"])),
#                   XM.norm=stringi::stri_sub(XM.norm, from=start-pos+1))
# )
# system.time(
#   bam.data <- base::rbind(bam.data.first,
#                           bam.data.second.leftmost,
#                           bam.data.second.rightmost,
#                           make.row.names=FALSE,
#                           stringsAsFactors=FALSE)
# )
# dim(bam.data)
# system.time(
#   base.counts <- sapply(XM.bases, stringi::stri_count_fixed, str=bam.data$XM.norm)
# )
# bam.data <- base::cbind(bam.data, base.counts)

# data.table is little faster
# system.time(
#     bam.data.merged <- data.table::setDT(bam.data)[,list(
#       start=min(pos, mpos),
#       XM=rcpp_merge_ends(flag, pos, isize, XM.norm)
#     ), by='qname,rname,XG']
# )

# # non-vectorised CX report
# system.time(
#   cx.report.nonvect <- bam.data.merged %>% #head(tail(bam.data.merged[order(bam.data.merged$start),], 20)) # bam.data.merged[1:1e4,] 
#     dplyr::mutate(rname=as.numeric(rname),
#                   strand=as.numeric(strand)) %>%
#     dplyr::select(rname, strand, start, XM, pass) %>%
#     apply(., 1, function (x) {
#       matrix(
#         rcpp_parse_xm(as.numeric(x[1]), as.numeric(x[2]), as.numeric(x[3]), x[4], as.logical(x[5]), ctx_str),
#         ncol=6, byrow=TRUE, dimnames=list(NULL, c("rname","strand","pos","ctx","meth","unmeth"))
#       )
#     }) %>%
#     do.call("rbind",.) %>%
#     as.data.frame() %>%
#     dplyr::group_by(rname,strand,pos,ctx) %>%
#     dplyr::summarise(meth=sum(meth), unmeth=sum(unmeth), .groups="drop") %>%
#     dplyr::mutate(rname=levels(bam.data.merged$rname)[rname],
#                   strand=levels(bam.data.merged$strand)[strand],
#                   ctx=intToUtf8(ctx, multiple=TRUE))
# )

# # fast merge reads, non-vectorised
# # Rcpp::sourceCpp("/Home/siv22/oni062/work/packages/epialleleR/src/rcpp_merge_ends.cpp")
# system.time(
#   if (alignment=="PE") {
#     bam.data.merged.nonvector <- bam.data %>%
#       dplyr::mutate(start=base::pmin.int(pos, mpos)) %>%
#       dplyr::group_by(qname, rname, strand=XG, start) %>%
#       dplyr::summarise(XM=rcpp_merge_ends(flag, pos, isize, XM.norm),
#                        .groups="drop")
#   } else {
#     bam.data.merged <- bam.data %>%
#       dplyr::select(qname, rname, strand=XG, start=pos, XM=XM.norm)
#   }
# )

# 
# # filtering before fast vectorised was introduced
# # gives NA when no out-of-context
# XM.bases       <- c("u","U","h","H","x","X","z","Z")
# system.time({
#   base.counts <- base::sapply(XM.bases, stringi::stri_count_fixed, str=bam.data.merged$XM)
#   class.counts <- as_tibble(sapply(XM.classes, function (cl) {matrixStats::rowSums2(base.counts[,cl,drop=FALSE])}))
#   bam.data.merged$pass <-
#     class.counts$ctx.meth+class.counts$ctx.unmeth >= min.context.sites &
#     class.counts$ctx.meth/(class.counts$ctx.meth+class.counts$ctx.unmeth) >= min.context.beta
#   if (filter.context!="any" & max.outofcontext.beta<1) {
#     bam.data.merged$pass <- bam.data.merged$pass &
#       class.counts$ooctx.meth/(class.counts$ooctx.meth+class.counts$ooctx.unmeth) <= max.outofcontext.beta
#   }
# })

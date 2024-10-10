#' plotPatterns
#'
#' @description
#' This convenience function plots methylation patterns (epialleles) previously
#' extracted by \code{\link[epialleleR]{extractPatterns}}.
#'
#' @details
#' As the number of methylation patterns can be quite large, by default, the
#' function plots \strong{the most abundant unique patterns} only. The complete
#' logic is as follows:
#' \itemize{
#'   \item from the input methylation patterns, all unique patterns are
#'   extracted and counted
#'   \item unique patterns are split in bins by their average beta value
#'   \item most abundant unique methylation patterns from each bin are plotted
#'   and silently returned
#' }
#' On the resulting plot, each cytosine is shown as a circle,
#' where the size of that circle represents cytosine context
#' and the fill encodes methylation status. If available, highlighted bases are
#' shown as labels of different colours.
#' 
#' @param patterns output of \code{\link[epialleleR]{preprocessBam}} function
#' (methylation patterns as a \code{\link[data.table]{data.table}} object).
#' @param order.by string defining order of patterns on the plot (default order
#' by: "beta").
#' @param beta.range numeric vector of length 2 for the range of average
#' pattern beta values represented on the plot (default: [0;1]).
#' @param bin.context string defining cytosine methylation context used
#' to calculate average beta value of a pattern that is further used to assign
#' patterns to bins:
#' \itemize{
#'   \item "CG" (the default) -- CpG cytosines (called as zZ)
#'   \item "CHG" -- CHG cytosines (xX)
#'   \item "CHH" -- CHH cytosines (hH)
#'   \item "CxG" -- CG and CHG cytosines (zZxX)
#'   \item "CX" -- all cytosines
#' }
#' @param nbins a single integer defining the number of bins (i.e., intervals
#' within `beta.range`). Default: 10.
#' @param npatterns.per.bin integer vector for the number of the most abundant
#' patterns selected from each bin (default: 2). When of length 1, the same
#' number of patterns will be taken. When of length `nbins`, allows to
#' fine-tune the number of selected patterns from each bin. Setting to `Inf`
#' effectively results in plotting all patterns.
#' @param plot.context string defining methylation context of cytosines
#' included in the plot (default: "CG"; for the range of available values, see
#' `bin.context` above).
#' @param genomic.scale string for the type of genomic position scale
#' of the plot: either "continuous" (the default) or "discrete".
#' @param breaks a vector of breaks for the genomic position scale of the plot.
#' If "auto" (the default), breaks for continuous scale are computed by the
#' default ggplot2 routines, while breaks for discrete scale are a subset of
#' plotted positions selected using \code{\link[base]{pretty}}.
#' Possible values: \code{\link[ggplot2:waiver]{ggplot2::waiver()}} for ggplot2
#' defaults, integer vector of breaks for continuous scale, or character vector
#' of breaks for discrete scale.
#' @param marginal string for the type of marginal plot: either "density"
#' (probability density of average beta values of all patterns; the default)
#' or "count" (counts of plotted patterns). "none" is not implemented yet; 
#' create an issue if interested.
#' @param marginal.position string for the position of marginal plot: either
#' "left" (the default) or "right" (not implemented yet; create an issue if
#' interested).
#' @param marginal.transform string for the transformation of marginal scale
#' (default: "identity"). Check
#' \code{\link[ggplot2:scale_x_continuous]{ggplot2::scale_x_continuous()}} for
#' more details.
#' @param marginal.limits limits of marginal scale (default: NULL). Check
#' \code{\link[ggplot2:scale_x_continuous]{ggplot2::scale_x_continuous()}} for
#' more details.
#' @param marginal.size numeric in range (0;1) for the relative width of the
#' marginal plot (default: 0.25).
#' @param ... additional arguments passed to
#' \code{\link[stats:density]{stats::density()}} call used in marginal density
#' plot. Possible value: \code{adjust=0.25}.
#' @param tag string for optional tagging of patterns with their count
#' ("count"), average beta value ("beta"), or pattern ID ("pattern").
#' Default: "none". 
#' @param tag.size numeric for the font size of the tag text
#' (in millimetres; default: 2.5).
#' @param tag.colour string for the colour of of the tag text.
#' Default: "#87654c".
#' @param tag.fill string for the colour of of the tag background.
#' Default: "lemonchiffon".
#' @param title the title of the plot. When `TRUE` (the default), a genomic
#' region from which patterns were extracted. Other possible values: anything
#' that can be converted to string, or `NULL` for no title.
#' @param subtitle the subtitle of the plot. When `TRUE` (the default), a
#' number of patterns plotted. Other possible values: anything
#' that can be converted to string, or `NULL` for no subtitle.
#' @param context.size a numeric vector with sizes of circles representing
#' cytosines within each of three contexts: CHH, CHG, and CG
#' (default: c(1, 2, 3)).
#' @param base.size numeric for the font size of the text for highlighted bases
#' (in millimetres; default: 3).
#' @param methylation.fill a vector of length 2 for colours representing
#' unmethylated and methylated cytosines, respectively. These colours are also
#' mapped to the lowest (0) and highest (1) possible beta values
#' to represent average beta values of methylation patterns and create a
#' gradient fill of a marginal density plot. Default: c("grey97", "grey10").
#' @param plot boolean. If `TRUE` (the default), patterns are plotted, and the
#' selected ones are silently returned as a \code{\link[data.table]{data.table}}
#' object. If `FALSE`, the \code{\link[gtable:gtable]{grob table}} object
#' is returned instead.
#' @param verbose boolean to report basic info on input and output.
#' @return the plot and (silently) the \code{\link[data.table]{data.table}}
#' object containing plotted methylation patterns (if `plot==TRUE`),
#' or \code{\link[gtable:gtable]{grob table}} object (if `plot==FALSE`).
#' @seealso \code{\link{extractPatterns}} for extracting methylation patterns,
#' \code{\link{preprocessBam}} for preloading BAM data,
#' \code{\link{generateCytosineReport}} for methylation statistics at the level
#' of individual cytosines, \code{\link{generateBedReport}} for genomic
#' region-based statistics, \code{\link{generateVcfReport}} for evaluating
#' epiallele-SNV associations, \code{\link{generateBedEcdf}} for analysing the
#' distribution of per-read beta values, and `epialleleR` vignettes for the
#' description of usage and sample data.
#' @examples
#'   # amplicon data
#'   amplicon.bam <- system.file("extdata", "amplicon010meth.bam",
#'                               package="epialleleR")
#'   custom.range <- as("chr17:43124861-43125150", "GRanges")
#'   
#'   # let's get our patterns
#'   patterns <- extractPatterns(bam=amplicon.bam, bed=custom.range)
#'   
#'   # default plot + silently returned plotted patterns
#'   selected.patterns <- plotPatterns(patterns)
#'   
#'   # all unique patterns with their counts as a margin, categorical positions,
#'   # tagged with pattern IDs, returned as a `gtable` object
#'   tbl <- plotPatterns(patterns, npatterns.per.bin=Inf, marginal="count",
#'                       genomic.scale="discrete", tag="pattern", plot=FALSE)
#'   
#'   # which can be plotted later
#'   grid::grid.newpage()
#'   grid::grid.draw(tbl)
#'   
#' @export
plotPatterns <- function (patterns, order.by=c("beta", "count"),
                          beta.range=c(0, 1), bin.context=c("CG", "CHG", "CHH", "CxG", "CX"), nbins=10, npatterns.per.bin=2,
                          plot.context=c("CG", "CHG", "CHH", "CxG", "CX"),
                          genomic.scale=c("continuous", "discrete"), breaks="auto",
                          marginal=c("density", "count"), marginal.position=c("left", "right"),
                          marginal.transform=c("identity", "log10"), marginal.limits=NULL, marginal.size=0.25, ...,
                          tag=c("none", "count", "beta", "pattern"), tag.size=2.5, tag.colour="#87654c", tag.fill="lemonchiffon",
                          title=TRUE, subtitle=TRUE, context.size=c(1, 2, 3), base.size=3, methylation.fill=c("grey97", "grey10"),
                          plot=TRUE, verbose=TRUE) {
  order.by <- match.arg(order.by)
  genomic.scale <- match.arg(genomic.scale)
  marginal.position <- match.arg(marginal.position)
  bin.context <- match.arg(bin.context)
  plot.context <- match.arg(plot.context)
  marginal <- match.arg(marginal)
  marginal.transform <- match.arg(marginal.transform)
  tag <- match.arg(tag)
  npatterns.per.bin <- rep(npatterns.per.bin, length.out=nbins)
  context.size <- rep(context.size, length.out=3)
  
  if (!requireNamespace("ggplot2", quietly=TRUE))
    stop("ggplot2 is required for plotting. Please install")
  
  if (plot.context=="CxG") {
    plot.context <- c("CG", "CHG")
  } else if (plot.context=="CX")  {
    plot.context <- c("CG", "CHG", "CHH")
  }
  bins <- seq(from=beta.range[1], to=beta.range[2], length.out=nbins+1)
  
  # all bases
  base.positions <- grep("^[0-9]+$", colnames(patterns), value=TRUE)
  patterns.summary <- patterns[, .(count=.N), by=c("pattern", base.positions)]
  
  context.to.factors <- lapply(.context.to.bases[[bin.context]], function (subctx) {
    bases <- unlist(strsplit(subctx, ""))
    match(bases, levels(patterns[[base.positions[1]]]))
  })
  
  patterns.summary[, beta:=apply(patterns.summary[, lapply(.SD, as.integer), .SDcols=base.positions], MARGIN=1, function (x) {
    meth <- sum(x %in% context.to.factors$ctx.meth, na.rm=TRUE)
    unmeth <- sum(x %in% context.to.factors$ctx.unmeth, na.rm=TRUE)
    beta <- meth/(meth+unmeth)
    return(if (meth+unmeth>0) beta else 0)
  })]
  
  patterns.summary[beta>=beta.range[1] & beta<=beta.range[2], bin:=findInterval(beta, bins, all.inside=TRUE)]
  patterns.selected <- patterns.summary[!is.na(bin)][order(-count), lapply(.SD, head, n=npatterns.per.bin[bin]), by=bin]
  patterns.selected[order(get(order.by), beta, count, decreasing=TRUE), I:=.N-.I]
  
  if (verbose) {
    bin.intervals <- levels(cut(bins, bins, include.lowest=TRUE, right=FALSE))
    nselected.per.bin <- data.table::merge.data.table(patterns.selected[, .N, by=bin], data.table::data.table(bin=seq_len(nbins)), all=TRUE)$N
    nselected.per.bin[is.na(nselected.per.bin)] <- 0
    stats <- sprintf(
      "%i patterns supplied\n%i unique\n%i most frequent unique patterns were selected for plotting using %i beta value bins:\n%s\n%s",
      nrow(patterns), nrow(patterns.summary), nrow(patterns.selected), nbins, paste(bin.intervals, collapse=" "),
      do.call("sprintf", c(list(fmt=paste(sprintf("%%%is", nchar(bin.intervals)), collapse=" ")), nselected.per.bin) )
    )
    message(stats)
  }
  
  plot.data <- data.table::melt.data.table(
    patterns.selected, measure.vars=base.positions, variable.name="pos", value.name="code", variable.factor=FALSE
  )[!is.na(code)]
  plot.data[, `:=` (
    base=factor(code, levels=c("A", "C", "G", "N", "T")),
    meth=factor(!code %in% c("h", "x", "z")),
    cntx=factor(tolower(code), levels=c("h", "x", "z"), labels=c("CHH", "CHG", "CG"))
  )]
  
  plot.data[, pos:=as.integer(pos)]
  if (genomic.scale=="discrete") {
    plot.data <- plot.data[is.na(cntx) | cntx %in% plot.context]
    plot.data[, pos:=factor(pos)]
  }
  if (identical(breaks, "auto")) {
    if (genomic.scale=="continuous") {
      breaks <- ggplot2::waiver()
    } else {
      breaks <- levels(plot.data$pos)[intersect(pretty(seq_along(levels(plot.data$pos))), seq_along(levels(plot.data$pos)))]
    }
  }
  
  # get title from bed
  if (identical(title, TRUE)) {
    title <- attr(patterns, "bed")
  }
  
  # get some subtitle stats
  if (identical(subtitle, TRUE)) {
    if (nrow(patterns.selected)==nrow(patterns.summary)) {
      subtitle <- sprintf("all %i unique patterns", nrow(patterns.selected))
    } else {
      subtitle <- sprintf("%i of %i unique patterns", nrow(patterns.selected), nrow(patterns.summary))
    }
  }
  
  # # marginal position
  # if (marginal.position=="right") {
  #   # todo
  # }
  
  # some kind of workaround for empty page but only when plotting
  if (plot) grid::grid.newpage()
  
  main.plot <- ggplot2::ggplot(plot.data, ggplot2::aes(x=pos, y=factor(I), group=factor(I))) +
    ggplot2::geom_line() +
    ggplot2::geom_segment(data=plot.data[, .(pos=sort(pos)[1]), by=I], mapping=ggplot2::aes(xend=-Inf, yend=factor(I)), linewidth=0.5, colour="grey") +
    ggplot2::geom_label(data=plot.data[!is.na(base)], mapping=ggplot2::aes(label=base, colour=base), size=base.size) +
    ggplot2::geom_point(data=plot.data[!is.na(cntx) & cntx %in% plot.context], mapping=ggplot2::aes(size=cntx, fill=meth), shape=21, colour=methylation.fill[2]) +
    ggplot2::scale_size_manual(name="context", values=setNames(context.size, c("CHH", "CHG", "CG"))) +
    ggplot2::scale_fill_manual(name="methylated", values=methylation.fill, drop=FALSE) +
    ggplot2::scale_y_discrete(name=NULL, breaks=NULL, labels=NULL) +
    ggplot2::theme_light() +
    ggplot2::theme(plot.margin=grid::unit(c(5.5, 5.5, 5.5, 0), "points")) +
    do.call(what=sprintf("scale_x_%s", genomic.scale), args=list(name="genomic position", breaks=breaks), envir=asNamespace("ggplot2"))
  
  # add tags if requested
  if (tag!="none") {
    tag.data <- plot.data[, .(pos=max(as.integer(pos))+0.5), by=.(I, label=get(tag))]
    if (tag=="beta") tag.data[, label:=sprintf("%.2g", label)]
    scale.range <- plot.data[, .(from=min(as.integer(pos)), to=max(as.integer(pos)))]
    tag.nchar <- max(nchar(as.character(tag.data$label)))
    tag.expand <- as.numeric(grid::convertX(grid::unit(tag.size*tag.nchar/2, "strwidth", "A"), "npc")) * (scale.range$to - scale.range$from + 1) + 1
    main.plot <- main.plot +
      ggplot2::geom_label(data=tag.data, mapping=ggplot2::aes(x=pos, y=factor(I), label=label), hjust=0, size=tag.size, colour=tag.colour, fill=tag.fill, inherit.aes=FALSE, show.legend=FALSE) +
      ggplot2::expand_limits(x=scale.range$to+tag.expand) +
      ggplot2::guides(tag=ggplot2::guide_custom(grid::legendGrob(tag, pch=22, do.lines=FALSE, gp=grid::gpar(col=tag.colour, fill=tag.fill, cex=3/3.88, lwd=0.75)) , title="tag"))
  }
  
  main.grob <- ggplot2::ggplotGrob(main.plot)
  
  if (marginal.transform=="identity") {
    trans <- "reverse"
    segment.x <- 0
    segment.xend <- -Inf
  } else {
    trans <- c("log10", "reverse")
    segment.x <- 1
    segment.xend <- 0
  }
  
  if (marginal=="count") { 
    side.plot <- ggplot2::ggplot(plot.data[, .(count=unique(count), beta=unique(beta)), by=I], ggplot2::aes(xmin=0, ymin=I-0.4, xmax=count, ymax=I+0.4, fill=beta)) +
      ggplot2::geom_rect(colour=methylation.fill[2]) +
      ggplot2::scale_x_continuous(transform=trans, limits=marginal.limits, minor_breaks=NULL, name="count") + # oob=scales::squish, 
      ggplot2::scale_y_continuous(name="patterns", breaks=NULL, minor_breaks=NULL, expand=ggplot2::expansion(0, 0.2)) + # , expand=ggplot2::expansion(0, max(plot.data$I)*0.03)) + # 
      ggplot2::scale_fill_gradient(low=methylation.fill[1], high=methylation.fill[2], limits=c(0, 1), guide="none") +
      ggplot2::theme_light() +
      ggplot2::theme(plot.margin=grid::unit(c(5.5, 0, 5.5, 5.5), "points"))
    marg.grob <- ggplot2::ggplotGrob(side.plot + ggplot2::ggtitle(title, subtitle=subtitle))
    marg.grob$widths[[marg.grob$layout[which(marg.grob$layout$name=="panel"), "l"]]] <- grid::unit(marginal.size/(1-marginal.size), "null")
  } else if (marginal=="density") {
    beta.dens <- data.table::as.data.table(
      stats::density(patterns.summary$beta, weights=patterns.summary$count/sum(patterns.summary$count),
                     from=min(patterns.summary$beta), to=max(patterns.summary$beta), n=2048, warnWbw=FALSE, ...)[c("x","y")]
    )
    plotted.min <- min(beta.dens[y>segment.xend]$y)
    beta.dens[y<plotted.min, y:=plotted.min]
    side.plot <- ggplot2::ggplot(beta.dens, ggplot2::aes(x=x, xend=x, y=y, yend=plotted.min, colour=x)) +
      ggplot2::geom_segment() + 
      ggplot2::geom_line(colour="black", linewidth=0.25) + 
      ggplot2::scale_colour_gradient(low=methylation.fill[1], high=methylation.fill[2], limits=c(0, 1), guide="none") +
      ggplot2::scale_y_continuous(transform=trans, limits=marginal.limits, minor_breaks=NULL, name="density") + # oob=scales::squish, 
      ggplot2::scale_x_continuous(name="average beta", breaks=bins, labels=format(bins, digits=2), minor_breaks=NULL, limits=beta.range, expand=ggplot2::expansion(0.05, 0)) +
      ggplot2::coord_flip() +
      ggplot2::geom_vline(mapping=ggplot2::aes(xintercept=x), data=data.table::data.table(x=bins), colour="white") +
      ggplot2::theme_light() +
      ggplot2::geom_segment(data=plot.data[, .(x=unique(beta))], mapping=ggplot2::aes(x=x, xend=x, y=min(beta.dens[y!=0]$y), yend=segment.xend), inherit.aes=FALSE, linewidth=0.5, colour="grey") +
      ggplot2::theme(plot.margin=grid::unit(c(5.5, 0, 5.5, 5.5), "points"))
    
    # corresponding lines
    max.one <- function (i, ...) { max(i, 1, ...) }
    corr.y <- function (i, b, mult=0.05) {
      b.scaled <- (b - beta.range[1]) / (beta.range[2] - beta.range[1])
      max.one(i) * (b.scaled + mult) / (1 + 2 * mult)
    }
    corr.ymax <- function (i, add=0.6) { (i+add) * max.one(i) / (max(i) + 2*add) }
    corr.plot <- ggplot2::ggplot(unique(plot.data[, .(I, beta)]), ggplot2::aes(x=1, xend=2, y=corr.y(I, beta), yend=corr.ymax(I))) +
      ggplot2::geom_segment(linewidth=0.5, colour="grey") +
      ggplot2::theme_void() +
      ggplot2::scale_x_continuous(expand=ggplot2::expansion(0, 0)) +
      ggplot2::scale_y_continuous(limits=c(0, max.one(plot.data$I)), expand=ggplot2::expansion(0, 0)) +
      ggplot2::theme(plot.margin=grid::unit(c(5.5, 0, 5.5, 0), "points"))

    side.grob <- ggplot2::ggplotGrob(side.plot + ggplot2::ggtitle(title, subtitle=subtitle))
    corr.grob <- ggplot2::ggplotGrob(corr.plot)
    side.grob$widths[[side.grob$layout[which(side.grob$layout$name=="panel"), "l"]]] <- grid::unit(marginal.size/(1-marginal.size), "null")
    corr.grob$widths[[corr.grob$layout[which(corr.grob$layout$name=="panel"), "l"]]] <- grid::unit(16, "points") # 12
    marg.grob <- cbind(side.grob, corr.grob)
  }
  
  comb.grob <- cbind(marg.grob, main.grob)
  for (i in which(comb.grob$layout$name=="background")) comb.grob$layout[i, c(1:4)] <- c(1, 1, dim(comb.grob))
  
  if (plot) {
    grid::grid.draw(comb.grob)
    return(invisible(patterns.selected))
  } else {
    return(comb.grob)
  }
}

# NB:
#   WORK IN PROGRESS
#
# TODO:
#   [x] add bed to attributes
#   [x] strand info is ignored in all methods - be clear on that
#   [ ] marginal=c("density", "count", "none")
#   [x] tag=c("none", "count", "beta", "pattern")
#   [x] tag guide
#   [x] proper expand for all plots
#   [x] scale labels too
#   [?] fill labels? N==NA? That would require that methylation aesthetics to be shape not colour
#   [x] import only required from ggplot? or fall back on its absence?
#   [x] make verbose work
#   [ ] make "right" work?




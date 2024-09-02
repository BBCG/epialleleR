# TODO:
#   [ ] add bed to attributes
#   [ ] marginal=c("density", "count", "none")
#   [ ] labels=c("none", "pattern", "count")
#   [ ] same expand for all plots, e.g., c(0, 0.5)
#   [ ] scale labels together with context?
#   [ ] fill labels? N==NA?
#   [ ] import only required from ggplot? or fall back on its absence?



plotPatterns <- function (patterns, order.by=c("beta", "count"),
                          beta.range=c(0, 1), bin.context=c("CG", "CHG", "CHH", "CxG", "CX"), nbins=10, max.patterns=20,
                          plot.context=c("CG", "CHG", "CHH", "CxG", "CX"),
                          genomic.scale=c("continuous", "discrete"), breaks="auto",
                          marginal=c("density", "count"), marginal.position=c("left", "right"), marginal.transform=c("identity", "log10"), marginal.limits=NULL, marginal.size=0.25, ...,
                          title=TRUE, subtitle=TRUE, context.size=1, colors=c("grey97", "grey10"),
                          plot=TRUE, verbose=TRUE) {
  order.by <- match.arg(order.by)
  genomic.scale <- match.arg(genomic.scale)
  marginal.position <- match.arg(marginal.position)
  bin.context <- match.arg(bin.context)
  plot.context <- match.arg(plot.context)
  marginal <- match.arg(marginal)
  marginal.transform <- match.arg(marginal.transform)
  
  # plot=TRUE; verbose=TRUE; order.by="beta"; genomic.scale="continuous"; breaks="auto"; beta.range=c(0, 1); bin.context="CG"; plot.context="CG"; nbins=10; max.patterns=20; title=TRUE; subtitle=TRUE; context.size=1; verbose=TRUE; marginal="density"; marginal.position="left"; marginal.transform="identity"; marginal.limits=NULL; colors=c("grey97", "grey10"); 
  # marginal="count"; genomic.scale="discrete"; marginal.transform="log10"
  require(ggplot2, quietly=TRUE)
  
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
  patterns.selected <- patterns.summary[!is.na(bin)][order(-count), lapply(.SD, head, n=max.patterns/nbins), by=bin]
  patterns.selected[order(get(order.by), beta, count, decreasing=TRUE), I:=.N-.I]
  
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
      breaks <- waiver()
    } else {
      breaks <- levels(plot.data$pos)[intersect(pretty(seq_along(levels(plot.data$pos))), seq_along(levels(plot.data$pos)))]
    }
  }
  
  # get title from bed
  if (title==TRUE) {
    title <- attr(patterns, "bed")
  }
  
  # get some subtitle stats
  if (subtitle==TRUE) {
    subtitle <- sprintf("%i of %i unique patterns", nrow(patterns.selected), nrow(patterns.summary))
  }
  
  # marginal position
  if (marginal.position=="right") {
    # todo
  }
  
  # need some kind of workaround for empty page
  
  
  main.plot <- ggplot(plot.data, aes(x=pos, y=I, group=I)) +
    geom_line() +
    geom_segment(data=plot.data[, .(pos=sort(pos)[1]), by=I], mapping=aes(xend=-Inf, yend=I), linewidth=0.5, colour="grey") +
    geom_label(data=plot.data[!is.na(base)], mapping=aes(label=base, color=base)) +
    geom_point(data=plot.data[!is.na(cntx) & cntx %in% plot.context], mapping=aes(size=cntx, fill=meth), shape=21, colour=colors[2]) +
    scale_size_manual(name="context", values=c("CHH"=1, "CHG"=2, "CG"=3)*context.size) +
    scale_fill_manual(name="methylated", values=colors) +
    scale_y_continuous(name=NULL, breaks=NULL, labels=NULL) +
    theme_light() +
    theme(plot.margin=grid::unit(c(5.5, 5.5, 5.5, 0), "points")) +
    do.call(sprintf("scale_x_%s", genomic.scale), list(name="genomic position", breaks=breaks))
  main.grob <- ggplotGrob(main.plot)
  
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
    # side.plot <- ggplot(plot.data[, .(count=unique(count)), by=I], aes(x=segment.x, y=I, xend=count, yend=I)) +
    #   geom_segment(color="lightgrey", linewidth=5) + #, ...) +
    #   scale_x_continuous(transform=trans, minor_breaks=NULL, name="count") +
    #   scale_y_continuous(name=NULL, breaks=NULL, minor_breaks=NULL) +
    #   theme_light() +
    #   geom_segment(data=plot.data[, .(y=unique(I))], mapping=aes(x=segment.x, xend=segment.xend, y=y, yend=y), inherit.aes=FALSE, linewidth=0.5, colour="grey") +
    #   theme(plot.margin=grid::unit(c(5.5, 0, 5.5, 5.5), "points"))
    # side.plot <- ggplot(plot.data[, .(count=unique(count), beta=unique(beta)), by=I], aes(xmin=segment.x, ymin=I-0.3, xmax=count, ymax=I+0.3, alpha=beta*alpha.scale+alpha.min)) +
    #   geom_rect(fill=meth.color, color=meth.color) + #, ...) +
    #   scale_x_continuous(transform=trans, minor_breaks=NULL, name="count") +
    #   scale_y_continuous(name=NULL, breaks=NULL, minor_breaks=NULL, expand=expansion(0,0.1)) +
    #   scale_alpha_identity(guide="none") +
    #   theme_light() +
    #   theme(plot.margin=grid::unit(c(5.5, 0, 5.5, 5.5), "points"))
    side.plot <- ggplot(plot.data[, .(count=unique(count), beta=unique(beta)), by=I], aes(xmin=0, ymin=I-0.3, xmax=count, ymax=I+0.3, fill=beta)) +
      geom_rect(color=colors[2]) + #, ...) +
      scale_x_continuous(transform=trans, minor_breaks=NULL, name="count") +
      scale_y_continuous(name="patterns", breaks=NULL, minor_breaks=NULL, expand=expansion(0, max(plot.data$I)*0.03)) + # , expand=expansion(0, 0) # 0.2/max(plot.data$I)
      scale_fill_gradient(low=colors[1], high=colors[2], limits=c(0, 1), guide="none") +
      theme_light() +
      theme(plot.margin=grid::unit(c(5.5, 0, 5.5, 5.5), "points"))
    marg.grob <- ggplotGrob(side.plot + ggtitle(title, subtitle=subtitle))
    marg.grob$widths[[marg.grob$layout[which(marg.grob$layout$name=="panel"), "l"]]] <- grid::unit(marginal.size/(1-marginal.size), "null")
  } else { # marginal=="density #
    # # plain density
    # side.plot <- ggplot(plot.data, aes(y=beta, weight=count)) +
    #   geom_density(fill="lightgrey") + #, ...) +
    #   scale_x_continuous(transform=trans, minor_breaks=NULL, name="density") +
    #   scale_y_continuous(breaks=bins, labels=format(bins, digits=2), minor_breaks=NULL, limits=beta.range) +
    #   geom_hline(mapping=aes(yintercept=y), data=data.table::data.table(y=bins), color="white") +
    #   theme_light() +
    #   geom_segment(data=plot.data[, .(y=unique(beta))], mapping=aes(x=segment.x, xend=segment.xend, y=y, yend=y), inherit.aes=FALSE, linewidth=0.5, colour="grey") +
    #   theme(plot.margin=grid::unit(c(5.5, 0, 5.5, 5.5), "points"))
    
    # # gradient density
    beta.dens <- data.table::as.data.table(
      density(patterns.summary$beta, weights=patterns.summary$count/sum(patterns.summary$count),
              from=min(patterns.summary$beta), to=max(patterns.summary$beta), n=2048, warnWbw=FALSE, ...)[c("x","y")]
    )
    plotted.min <- min(beta.dens[y>segment.xend]$y)
    beta.dens[y<plotted.min, y:=plotted.min] # beta.dens[y<plotted.min, y:=NA] # 
    side.plot <- ggplot(beta.dens, aes(x=x, xend=x, y=y, yend=plotted.min, color=x)) +
      geom_segment() + 
      geom_line(color="black", linewidth=0.25) + 
      scale_color_gradient(low=colors[1], high=colors[2], limits=c(0, 1), guide="none") +
      scale_y_continuous(transform=trans, minor_breaks=NULL, name="density") +
      scale_x_continuous(name="average beta", breaks=bins, labels=format(bins, digits=2), minor_breaks=NULL, limits=beta.range) +
      coord_flip() +
      geom_vline(mapping=aes(xintercept=x), data=data.table::data.table(x=bins), color="white") +
      theme_light() +
      geom_segment(data=plot.data[, .(x=unique(beta))], mapping=aes(x=x, xend=x, y=min(beta.dens[y!=0]$y), yend=segment.xend), inherit.aes=FALSE, linewidth=0.5, colour="grey") +
      theme(plot.margin=grid::unit(c(5.5, 0, 5.5, 5.5), "points"))
    # # marginal=="histogram"
    #   side.plot <- ggplot(plot.data, aes(y=beta, weight=count)) +
    #     geom_histogram(fill="lightgrey", boundary=0) + #, ...) +
    #     scale_x_continuous(transform=trans, minor_breaks=NULL, name="histogram") + 
    #     scale_y_continuous(breaks=bins, labels=format(bins, digits=2), minor_breaks=NULL, limits=beta.range) +
    #     geom_hline(mapping=aes(yintercept=y), data=data.table::data.table(y=bins), color="white") +
    #     theme_light() +
    #     geom_segment(data=plot.data[, .(y=unique(beta))], mapping=aes(x=segment.x, xend=segment.xend, y=y, yend=y), inherit.aes=FALSE, linewidth=0.5, colour="grey") +
    #     theme(plot.margin=grid::unit(c(5.5, 0, 5.5, 5.5), "points"))
    
    # corresponding lines
    corr.plot <- ggplot(unique(plot.data[, .(I, beta)]), aes(x=1, xend=2, y=max(I)*(beta-beta.range[1])/(beta.range[2]-beta.range[1]), yend=I)) +
      geom_segment(linewidth=0.5, colour="grey") +
      theme_void() +
      scale_x_continuous(expand=expansion(0,0)) +
      theme(plot.margin=grid::unit(c(5.5, 0, 5.5, 0), "points"))
    
    side.grob <- ggplotGrob(side.plot + ggtitle(title, subtitle=subtitle))
    corr.grob <- ggplotGrob(corr.plot)
    side.grob$widths[[side.grob$layout[which(side.grob$layout$name=="panel"), "l"]]] <- grid::unit(marginal.size/(1-marginal.size), "null")
    corr.grob$widths[[corr.grob$layout[which(corr.grob$layout$name=="panel"), "l"]]] <- grid::unit(12, "points")
    marg.grob <- cbind(side.grob, corr.grob)
  }
  
  comb.grob <- cbind(marg.grob, main.grob)
  for (i in which(comb.grob$layout$name=="background")) comb.grob$layout[i, c(1:4)] <- c(1, 1, dim(comb.grob))
  
  if (plot) {
    plot(comb.grob)
    return(invisible(patterns.selected))
  } else {
    return(comb.grob)
  }
}




#########################################################################/**
# @RdocFunction plot
#
# @alias plot,QDNAseqSignals,missing-method
#
# @title "Plot copy number profile"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{x}{...}
#   \item{y}{...}
#   \item{...}{...}
# }
#
# \value{
#   Returns a named @list containing elements ...
# }
#
# @author "IS"
#
# \seealso{
#   Internally, ...
# }
#
#*/#########################################################################
setMethod("plot", signature(x="QDNAseqSignals", y="missing"),
  function (x, y, main=NULL, includeReadCounts=TRUE,
  logTransform=TRUE, sdFUN=sdDiff,
  delcol="darkred", losscol="red", gaincol="blue", ampcol="darkblue",
  pointcol="black", segcol="chocolate", misscol=NA,
  ylab=NULL, ylim=NULL, yaxp=NULL, ... ) {
  if (class(x) %in% c("QDNAseqCopyNumbers", "QDNAseqReadCounts")) {
    condition <- binsToUse(x)
  } else {
    condition <- rep(TRUE, times=nrow(x))
  }
  baseLine <- NA
  if ("calls" %in% assayDataElementNames(x)) {
    if (is.null(ylim))
      if (logTransform) {
        ylim <- c(-5, 5)
      } else {
        ylim <- c(-2, 4)
      }
  }
  if ("copynumber" %in% assayDataElementNames(x)) {
    copynumber <- assayDataElement(x, "copynumber")[condition, , drop=FALSE]
    if (is.null(ylab))
      ylab <- ifelse(logTransform, expression(log[2]~ratio), "ratio")
    if (is.null(ylim))
      if (logTransform) {
        ylim <- c(-3, 5)
      } else {
        ylim <- c(0, 4)
      }
    if (is.null(yaxp))
      yaxp <- c(ylim[1], ylim[2], ylim[2]-ylim[1])
    baseLine <- ifelse(logTransform, 0, 1)
  } else {
    copynumber <- assayDataElement(x, "counts")[condition, , drop=FALSE]
    if (is.null(ylab))
      ylab <- ifelse(logTransform, expression(log[2]~read~count),
        "read count")
    if (is.null(ylim))
      if (logTransform) {
        ylim <- c(0, max(log2adhoc(copynumber)))
      } else {
        ylim <- range(copynumber)
      }
  }
  sdFUN <- match.fun(sdFUN)
  noise <- apply(scale(copynumber, center=FALSE,
    scale=apply(copynumber, 2, mean, na.rm=TRUE)), 2, sdFUN, na.rm=TRUE)
  if (logTransform)
    copynumber <- log2adhoc(copynumber)
  if (is.null(main))
    main <- sampleNames(x)
  if (includeReadCounts && "reads" %in% names(pData(x)))
    main <- paste(main, " (", format(x$reads, trim=TRUE, big.mark=","),
      " reads)", sep="")
  if (length(ylab) == 1)
    ylab <- rep(ylab, times=ncol(x))
  all.chrom <- chromosomes(x)
  if (class(x) %in% c("cghRaw", "cghSeg", "cghCall")) {
    chrom.lengths <- CGHbase:::.getChromosomeLengths("GRCh37")
  } else {
    all.chrom.lengths <- aggregate(bpend(x),
      by=list(chromosome=all.chrom), FUN=max)
    chrom.lengths <- all.chrom.lengths$x
    names(chrom.lengths) <- all.chrom.lengths$chromosome
  }
  chrom <- all.chrom[condition]
  uni.chrom <- unique(chrom)
  chrom.lengths <- chrom.lengths[as.character(uni.chrom)]
  pos <- as.numeric(bpstart(x)[condition])
  pos2 <- as.numeric(bpend(x)[condition])
  chrom.ends <- integer()
  cumul <- 0
  for (i in uni.chrom) {
    pos[chrom > i] <- pos[chrom > i] + chrom.lengths[as.character(i)]
    pos2[chrom > i] <- pos2[chrom > i] + chrom.lengths[as.character(i)]
    cumul <- cumul + chrom.lengths[as.character(i)]
    chrom.ends <- c(chrom.ends, cumul)
  }
  names(chrom.ends) <- names(chrom.lengths)
  for (i in seq_len(ncol(x))) {
    vmsg("Plotting sample ", main[i])
    cn <- copynumber[, i]
    if ("segmented" %in% assayDataElementNames(x)) {
      segmented <- assayDataElement(x, "segmented")[condition, i]
      if (logTransform)
        segmented <- log2adhoc(segmented)
      segment <- CGHbase:::.makeSegments(segmented, chrom)
    }
    if ("calls" %in% assayDataElementNames(x)) {
      losses <- probloss(x)[condition, i]
      gains <- probgain(x)[condition, i]
      if (!is.null(probdloss(x)))
        losses <- losses + probdloss(x)[condition, i]
      if (!is.null(probamp(x)))
        gains <- gains + probamp(x)[condition, i]
      par(mar=c(5, 4, 4, 4) + 0.2)
      plot(NA, main=main[i], xlab=NA, ylab=NA, las=1,
        xlim=c(0, max(pos2)), ylim=ylim, xaxs="i", xaxt="n",
        yaxp=c(ylim[1], ylim[2], ylim[2]-ylim[1]), yaxs="i")
      lim <- par("usr")
      lim[3:4] <- c(0, 1)
      par(usr=lim)
      dticks <- seq(0, 1, by=0.2)
      axis(4, at=dticks, labels=NA, srt=270, las=1, cex.axis=1,
        cex.lab=1, tck=-0.015)
      axis(4, at=dticks, labels=dticks, srt=270, las=1, cex.axis=1,
        cex.lab=1, line=-0.4, lwd=0)
      mtext("probability", side=4, line=2, cex=par("cex"))
      if (!is.na(misscol)) {
        rect(0, -1, max(pos2), 1, col=misscol, border=NA)
        rect(pos, -1, pos2, 1, col="white", border=NA)
      }
      rect(pos[segment[,2]], 0, pos2[segment[,3]], losses[segment[,2]],
        col=losscol, border=losscol)
      rect(pos[segment[,2]], 1, pos2[segment[,3]], 1-gains[segment[,2]],
        col=gaincol, border=gaincol)
      axis(3, at=pos[which(probamp(x)[condition,i] >= 0.5)],
        labels=FALSE, col=ampcol, col.axis="black", srt=270, las=1,
        cex.axis=1, cex.lab=1)
      axis(1, at=pos[which(probdloss(x)[condition,i] >= 0.5)],
        labels=FALSE, col=delcol, col.axis="black", srt=270, las=1,
        cex.axis=1, cex.lab=1)
      box()
      lim[3:4] <- ylim
      par(usr=lim)
      points(pos, cn, cex=.1, col=pointcol)
    } else {
      plot(pos, cn, cex=.1, col=pointcol, main=main[i],
        xlab=NA, ylab=NA, ylim=ylim, xaxt="n", xaxs="i", yaxs="i",
        yaxp=yaxp, tck=-0.015, las=1)
    }
    mtext(text="chromosomes", side=1, line=2, cex=par("cex"))
    mtext(text=ylab[i], side=2, line=2, cex=par("cex"))
    abline(h=baseLine)
    abline(v=chrom.ends[-length(chrom.ends)], lty="dashed")
    ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
    axis(side=1, at=ax, labels=NA, cex=.2, lwd=.5, las=1, cex.axis=1,
      cex.lab=1, tck=-0.015)
    axis(side=1, at=ax, labels=uni.chrom, cex=.2, lwd=0, las=1, cex.axis=1,
      cex.lab=1, tck=-0.015, line=-0.4)
    if ("segmented" %in% assayDataElementNames(x)) {
      for (jjj in seq_len(nrow(segment))) {
        segments(pos[segment[jjj,2]], segment[jjj,1], pos[segment[jjj,3]],
          segment[jjj,1], col=segcol, lwd=3)
      }
    }
    par(xpd=TRUE)
    amps <- cn
    amps[amps <= ylim[2]] <- NA_real_
    amps[!is.na(amps)] <- ylim[2] + 0.01 * (ylim[2]-ylim[1])
    dels <- cn
    dels[dels >= ylim[1]] <- NA_real_
    dels[!is.na(dels)] <- ylim[1] - 0.01 * (ylim[2]-ylim[1])
    points(pos, amps, pch=24, col=pointcol, bg=pointcol, cex=0.5)
    points(pos, dels, pch=25, col=pointcol, bg=pointcol, cex=0.5)
    if ("segmented" %in% assayDataElementNames(x)) {
      amps <- assayDataElement(x, "segmented")[condition, i]
      amps[amps <= ylim[2]] <- NA_real_
      amps[!is.na(amps)] <- ylim[2] + 0.01 * (ylim[2]-ylim[1])
      dels <- assayDataElement(x, "segmented")[condition, i]
      dels[dels >= ylim[1]] <- NA_real_
      dels[!is.na(dels)] <- ylim[1] - 0.01 * (ylim[2]-ylim[1])
      points(pos, amps, pch=24, col=segcol, bg=segcol, cex=0.5)
      points(pos, dels, pch=25, col=segcol, bg=segcol, cex=0.5)
    }
    par(xpd=FALSE)
    ### estimate for standard deviation
    if ("generic" %in% slotNames(sdFUN) &&
      length(grep("Diff$", sdFUN@generic)) == 1) {
      symbol <- quote(hat(sigma)[Delta])
    } else {
      symbol <- substitute(hat(sigma))
    }
    mtext(substitute(symbol==sd, list(symbol=symbol, sd=sprintf("%.3g",
      noise[i]))), side=3, line=0, adj=1, cex=par("cex"))
    ### number of data points
    str <- paste(round(sum(condition) / 1000), "k x ", sep="")
    probe <- median(bpend(x)-bpstart(x)+1)
    if (probe < 1000) {
      str <- paste(str, probe, " bp", sep="")
    } else {
      str <- paste(str, round(probe / 1000), " kbp", sep="")
    }
    if ("segmented" %in% assayDataElementNames(x))
      str <- paste(str, ", ", nrow(segment), " segments", sep="")
    mtext(str, side=3, line=0, adj=0, cex=par("cex"))
  }
})




#########################################################################/**
# @RdocFunction frequencyPlot
#
# @alias frequencyPlot,QDNAseqCopyNumbers,missing-method
#
# @title "Plot copy number aberration frequencies"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{x}{...}
#   \item{y}{...}
#   \item{...}{...}
# }
#
# \value{
#   Returns a named @list containing elements ...
# }
#
# @author "IS"
#
# \seealso{
#   Internally, ...
# }
#
#*/#########################################################################
setMethod("frequencyPlot", signature=c(x="QDNAseqCopyNumbers", y="missing"),
  function(x, y, main="Frequency Plot", losscol="red", gaincol="blue",
  misscol=NA, ... ) {
  condition <- binsToUse(x)
  all.chrom <- chromosomes(x)
  all.chrom.lengths <- aggregate(bpend(x),
    by=list(chromosome=all.chrom), FUN=max)
  chrom.lengths <- all.chrom.lengths$x
  names(chrom.lengths) <- all.chrom.lengths$chromosome
  chrom <- all.chrom[condition]
  uni.chrom <- unique(chrom)
  chrom.lengths <- chrom.lengths[as.character(uni.chrom)]
  pos <- as.numeric(bpstart(x)[condition])
  pos2 <- as.numeric(bpend(x)[condition])
  chrom.ends <- integer()
  cumul <- 0
  for (i in uni.chrom) {
    pos[chrom > i] <- pos[chrom > i] + chrom.lengths[as.character(i)]
    pos2[chrom > i] <- pos2[chrom > i] + chrom.lengths[as.character(i)]
    cumul <- cumul + chrom.lengths[as.character(i)]
    chrom.ends <- c(chrom.ends, cumul)
  }
  names(chrom.ends) <- names(chrom.lengths)
  calls <- calls(x)[condition, , drop=FALSE]
  loss.freq <- rowMeans(calls < 0)
  gain.freq <- rowMeans(calls > 0)
  plot(NA, main=main, xlab="chromosomes", ylab="frequency",
    xlim=c(0, max(pos2)), ylim=c(-1,1), xaxs="i", xaxt="n",
    yaxs="i", yaxt="n",...)
  if (!is.na(misscol)) {
    rect(0, -1, max(pos2), 1, col=misscol, border=NA)
    rect(pos, -1, pos2, 1, col="white", border=NA)
  }
  rect(pos, 0, pos2, gain.freq, col=gaincol, border=gaincol)
  rect(pos, 0, pos2, -loss.freq, col=losscol, border=losscol)
  box()
  abline(h=0)
  abline(v=chrom.ends[-length(chrom.ends)], lty="dashed")
  ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
  axis(side=1, at=ax, labels=uni.chrom, cex=.2, lwd=.5, las=1, cex.axis=1,
    cex.lab=1)
  axis(side=2, at=c(-1, -0.5, 0, 0.5, 1), labels=c("100 %", " 50 %", "0 %",
    "50 %", "100 %"), las=1)
  mtext("gains", side=2, line=3, at=0.5, cex=par("cex"))
  mtext("losses", side=2, line=3, at=-0.5, cex=par("cex"))
  ### number of data points
  str <- paste(round(sum(condition) / 1000), "k x ", sep="")
  probe <- median(bpend(x)-bpstart(x)+1)
  if (probe < 1000) {
    str <- paste(str, probe, " bp", sep="")
  } else {
    str <- paste(str, round(probe / 1000), " kbp", sep="")
  }
  mtext(str, side=3, line=0, adj=0, cex=par("cex"))
})




#########################################################################/**
# @RdocFunction isobarPlot
#
# @alias isobarPlot,QDNAseqReadCounts,missing-method
#
# @title "Plot median read counts as a function of GC content and mappability"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{x}{...}
#   \item{y}{...}
#   \item{...}{...}
# }
#
# \value{
#   Returns a named @list containing elements ...
# }
#
# @author "IS"
#
# \seealso{
#   Internally, ...
# }
#
#*/#########################################################################
setMethod("isobarPlot", signature=c(x="QDNAseqReadCounts", y="missing"),
  definition=function(x, y, main=NULL,
  what=c("read counts", "fit", "residuals"), adjustIncompletes=TRUE,
  ...) {
  what <- match.arg(what)
  if (what == "read counts") {
    signal <- assayDataElement(x, "counts")
    if (adjustIncompletes) {
      signal <- signal / fData(x)$bases * 100L
      signal[fData(x)$bases == 0] <- 0L
    }
    if (is.null(main))
      main <- paste(sampleNames(x), "median read counts")
  } else if (what == "fit") {
    signal <- assayDataElement(x, "fit")
    if (is.null(main))
      main <- paste(sampleNames(x), "loess fit")
  } else if (what == "residuals") {
    counts <- assayDataElement(x, "counts")
    if (adjustIncompletes) {
      counts <- counts / fData(x)$bases * 100L
      counts[fData(x)$bases == 0] <- 0L
    }
    fit <- assayDataElement(x, "fit")
    signal <- counts / fit
    signal[fit <= 0] <- 0
    if (is.null(main))
      main <- paste(sampleNames(x), "median loess residuals")
  }
  condition <- binsToUse(x)
  gc <- round(fData(x)$gc)
  mappability <- round(fData(x)$mappability)
  median.signal <- aggregate(signal[condition, ], by=list(gc=gc[condition],
    mappability=mappability[condition]), FUN=median)
  median.signal <- median.signal[!is.na(median.signal$gc), ]
  median.signal <- median.signal[!is.na(median.signal$mappability), ]
  rownames(median.signal) <- paste(median.signal$gc, "-",
    median.signal$mappability, sep="")
  xx <- min(median.signal$mappability):max(median.signal$mappability)
  yy <- min(median.signal$gc):max(median.signal$gc)
  m <- matrix(nrow=length(xx), ncol=length(yy), dimnames=list(xx, yy))

  for (i in seq_len(ncol(x))) {
    vmsg("Plotting sample ", main[i])
    for (j in 1:nrow(median.signal))
      m[as.character(median.signal[j, "mappability"]),
        as.character(median.signal[j, "gc"])] <- median.signal[j, i+2L]
    image(xx, yy, m, col=paste("#", c(sprintf("%02X", 0L:255L),
      rep("FF", 256L)), c(rep("FF", 256L), sprintf("%02X", 255L:0L)),
      sprintf("%02X", 255L), sep=""),
      main=main[i],
      xlab=NA, ylab=NA, xaxt="n", yaxt="n",
      zlim=c(-max(abs(range(m, finite=TRUE))),
      max(abs(range(m, finite=TRUE)))), ...)
    axis(side=1, tck=-.015, labels=NA)
    axis(side=2, tck=-.015, labels=NA)
    axis(side=1, lwd=0, line=-0.4)
    axis(side=2, lwd=0, line=-0.4)
    mtext(side=1, "mappability", line=2, cex=par("cex"))
    mtext(side=2, "GC content", line=2, cex=par("cex"))
    contour(xx, yy, m, nlevels=20L, zlim=c(-max(abs(range(m, finite=TRUE))),
      max(abs(range(m, finite=TRUE)))), add=TRUE)

    str <- paste(round(sum(condition) / 1000), "k x ", sep="")
    probe <- median(bpend(x)-bpstart(x)+1)
    if (probe < 1000) {
      str <- paste(str, probe, " bp", sep="")
    } else {
      str <- paste(str, round(probe / 1000), " kbp", sep="")
    }
    mtext(str, side=3, line=0, adj=0, cex=par("cex"))

    reads <- sum(assayDataElement(x, "counts")[condition, i], na.rm=TRUE)
    mtext(paste(format(reads, trim=TRUE, big.mark=","), " reads",
      sep=""), side=3, line=0, adj=1, cex=par("cex"))
  }
})




#########################################################################/**
# @RdocFunction noisePlot
#
# @alias noisePlot,QDNAseqReadCounts,missing-method
#
# @title "Plot noise as a function of sequence depth"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{x}{A @see "QDNAseqReadCounts" object.}
#   \item{y}{missing}
#   \item{...}{Further arguments to @see "graphics::plot" and
#     @see "graphics::text".}
# }
#
# @author "IS"
#
#*/#########################################################################
setMethod("noisePlot", signature=c(x="QDNAseqReadCounts", y="missing"),
  definition=function(x, y, main="Noise Plot", fit=NULL, ...) {
  condition <- binsToUse(x)
  totalReads <- x$reads
  counts <- assayDataElement(x, "counts")[condition, , drop=FALSE]
  usedReads <- apply(counts, 2, sum)
  reciprocalOfAverageUsedReadsPerBin <- 1/(usedReads/sum(condition))
  if (is.null(fit)) {
    if (! "fit" %in% assayDataElementNames(x))
      x <- estimateCorrection(x)
    fit <- assayDataElement(x, "fit")[condition, , drop=FALSE]
  }
  signal <- counts / fit
  signal[fit <= 0] <- 0
  signal <- scale(signal, center=FALSE,
    scale=apply(signal, 2, mean, na.rm=TRUE))
  noise <- apply(signal, 2, sdDiff, na.rm=TRUE)
  plot(reciprocalOfAverageUsedReadsPerBin, noise^2, main=main, cex=0.5,
    xlim=c(0, 1.1*max(reciprocalOfAverageUsedReadsPerBin)),
    ylim=c(0, 1.04*max(noise^2)),
    xlab=NA, ylab=NA, xaxs="i", yaxs="i", xaxt="n", yaxt="n", ...)
  if (ncol(x) > 1) {
    relationship <- lm(totalReads ~ usedReads)
  } else {
    relationship <- lm(totalReads ~ 0 + usedReads)
  }
  at <- axTicks(side=1)
  usedReadsAtTicks <- sum(condition)/at
  labels <- round(predict(relationship,
    newdata=data.frame(usedReads=usedReadsAtTicks)) / 1e6, digits=1)
  labels[1] <- NA
  axis(side=1, tck=-.015, at=at, labels=NA)
  axis(side=2, tck=-.015, labels=NA)
  axis(side=1, lwd=0, line=-0.4, at=at, labels=labels)
  axis(side=2, lwd=0, line=-0.4)
  mtext(side=1, "million reads", line=2, cex=par("cex"))
  mtext(side=2, expression(hat(sigma)[Delta]^2), line=2, las=1, cex=par("cex"))
  text(reciprocalOfAverageUsedReadsPerBin, noise^2,
    labels=sampleNames(x), pos=4, cex=0.5, ...)
  abline(0, 1)
})

# EOF

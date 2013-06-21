#########################################################################/**
# @RdocFunction plot
#
# @alias plot,QDNAseqReadCounts,missing-method
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
setMethod('plot', signature(x='QDNAseqReadCounts', y='missing'),
  function (x, y, main=NULL, includeReadCounts=TRUE,
    delcol='darkred', losscol='red', gaincol='blue', ampcol='darkblue',
    pointcol='black', segcol='chocolate', misscol=NA,
    ylab=expression(normalized~log[2]~read~count), ylim=NULL, ... ) {
    if (is.null(ylim)) {
      if ('calls' %in% assayDataElementNames(x)) {
        ylim <- c(-5, 5)
      } else {
        ylim <- c(-2, 5)
      }
    }
    if (is.null(main))
      main <- sampleNames(x)
   if (includeReadCounts && 'reads' %in% names(pData(x)))
      main <- paste(main, ' (', format(x$reads, trim=TRUE, big.mark=','),
        ' reads)', sep='')
    if ('filter' %in% colnames(fData(x))) {
      condition <- fData(x)$filter
    } else {
      conditon <- rep(TRUE, times=nrow(x))
    }
    all.chrom <- chromosomes(x)
    all.chrom.lengths <- aggregate(bpend(x),
      by=list(chromosome=all.chrom), max)
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
    copynumber <- copynumber(x)[condition, , drop=FALSE]
    for (i in 1:ncol(x)) {
      message('Plotting sample ', main[i])
      if ('segmented' %in% assayDataElementNames(x))
        segment <- CGHbase:::.makeSegments(segmented(x)
          [condition, i], chrom)
      if ('calls' %in% assayDataElementNames(x)) {
        losses <- probloss(x)[condition, i]
        gains <- probgain(x)[condition, i]
        if (!is.null(probdloss(x)))
          losses <- losses + probdloss(x)[condition, i]
        if (!is.null(probamp(x)))
          gains <- gains + probamp(x)[condition, i]
        par(mar=c(5, 4, 4, 4) + 0.2)
        plot(NA, main=main[i], xlab='chromosomes', ylab=ylab, las=1,
          xlim=c(0, max(pos2)), ylim=ylim, xaxs='i', xaxt='n',
          yaxp=c(ylim[1], ylim[2], ylim[2]-ylim[1]), yaxs='i')
        lim <- par('usr')
        lim[3:4] <- c(0, 1)
        par(usr=lim)
        dticks <- seq(0, 1, by=0.2)
        axis(4, at=dticks, labels=dticks, srt=270, las=1, cex.axis=1,
          cex.lab=1)
        mtext('probability', side=4, line=3, cex=par('cex'))
        if (!is.na(misscol)) {
          rect(0, -1, max(pos2), 1, col=misscol, border=NA)
          rect(pos, -1, pos2, 1, col='white', border=NA)
        }
        rect(pos[segment[,2]], 0, pos2[segment[,3]], losses[segment[,2]],
          col=losscol, border=losscol)
        rect(pos[segment[,2]], 1, pos2[segment[,3]], 1-gains[segment[,2]],
          col=gaincol, border=gaincol)
        axis(3, at=pos[which(probamp(x)[condition,i] >= 0.5)],
          labels=FALSE, col=ampcol, col.axis='black', srt=270, las=1,
          cex.axis=1, cex.lab=1)
        axis(1, at=pos[which(probdloss(x)[condition,i] >= 0.5)],
          labels=FALSE, col=delcol, col.axis='black', srt=270, las=1,
          cex.axis=1, cex.lab=1)
        box()
        lim[3:4] <- ylim
        par(usr=lim)
        points(pos, copynumber[, i], cex=.1, col=pointcol)
      } else {
        plot(pos, copynumber[, i], cex=.1, col=pointcol, main=main[i],
          xlab='chromosomes', ylab=ylab, ylim=ylim, xaxt='n', xaxs='i',
          yaxp=c(ylim[1], ylim[2], ylim[2]-ylim[1]), yaxs='i')
      }
      abline(h=0)
      abline(v=chrom.ends[-length(chrom.ends)], lty='dashed')
      ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
      axis(side=1, at=ax, labels=uni.chrom, cex=.2, lwd=.5, las=1, cex.axis=1,
        cex.lab=1)
      if ('segmented' %in% assayDataElementNames(x)) {
        for (jjj in (1:nrow(segment))) {
          segments(pos[segment[jjj,2]], segment[jjj,1], pos[segment[jjj,3]],
            segment[jjj,1], col=segcol, lwd=3)        
        }
      }
      amps <- copynumber[, i]
      amps[amps < ylim[2]] <- NA
      amps[!is.na(amps)] <- ylim[2] + 0.01 * (ylim[2]-ylim[1])
      dels <- copynumber[, i]
      dels[dels > ylim[1]] <- NA
      dels[!is.na(dels)] <- ylim[1] - 0.01 * (ylim[2]-ylim[1])
      par(xpd=TRUE)
      points(pos, amps, pch=24, col=pointcol, bg=pointcol, cex=0.5)
      points(pos, dels, pch=25, col=pointcol, bg=pointcol, cex=0.5)
      par(xpd=FALSE)
      ### MAD
      mtext(substitute(hat(sigma)[Delta]==sd, list(sd=sprintf('%.3g',
        madDiff(copynumber[, i], na.rm=TRUE)))), side=3, line=0,
        adj=1, cex=par('cex'))
      ### number of data points
      str <- paste(round(length(chrom) / 1000), 'k x ', sep='')
      probe <- median(bpend(x)-bpstart(x)+1)
      if (probe < 1000) {
        str <- paste(str, probe, ' bp', sep='')
      } else {
        str <- paste(str, round(probe / 1000), ' kbp', sep='')
      }
      if ('segmented' %in% assayDataElementNames(x))
        str <- paste(str, ', ', nrow(segment), ' segments', sep='')
      mtext(str, side=3, line=0, adj=0, cex=par('cex'))
    }
})




#########################################################################/**
# @RdocFunction frequencyPlot
#
# @alias frequencyPlot,QDNAseqReadCounts,missing-method
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
setMethod('frequencyPlot', signature=c(x='QDNAseqReadCounts', y='missing'),
  function(x, y, main='Frequency Plot', losscol='red', gaincol='blue',
  misscol=NA, ... ) {
  if ('filter' %in% colnames(fData(x))) {
    condition <- fData(x)$filter
  } else {
    conditon <- rep(TRUE, times=nrow(x))
  }
  all.chrom <- chromosomes(x)
  all.chrom.lengths <- aggregate(bpend(x),
    by=list(chromosome=all.chrom), max)
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
  plot(NA, main=main, xlab='chromosomes', ylab='frequency',
    xlim=c(0, max(pos2)), ylim=c(-1,1), xaxs='i', xaxt='n',
    yaxs='i', yaxt='n',...)
  if (!is.na(misscol)) {
    rect(0, -1, max(pos2), 1, col=misscol, border=NA)
    rect(pos, -1, pos2, 1, col='white', border=NA)
  }
  rect(pos, 0, pos2, gain.freq, col=gaincol, border=gaincol)
  rect(pos, 0, pos2, -loss.freq, col=losscol, border=losscol)
  box()
  abline(h=0)
  abline(v=chrom.ends[-length(chrom.ends)], lty='dashed')
  ax <- (chrom.ends + c(0, chrom.ends[-length(chrom.ends)])) / 2
  axis(side=1, at=ax, labels=uni.chrom, cex=.2, lwd=.5, las=1, cex.axis=1,
    cex.lab=1)
  axis(side=2, at=c(-1, -0.5, 0, 0.5, 1), labels=c('100 %', ' 50 %', '0 %',
    '50 %', '100 %'), las=1)
  mtext('gains', side=2, line=3, at=0.5)
  mtext('losses', side=2, line=3, at=-0.5)
  ### number of data points
  str <- paste(round(nrow(x) / 1000), 'k x ', sep='')
  probe <- median(bpend(x)-bpstart(x)+1)
  if (probe < 1000) {
    str <- paste(str, probe, ' bp', sep='')
  } else {
    str <- paste(str, round(probe / 1000), ' kbp', sep='')
  }
  mtext(str, side=3, line=0, adj=0)
})

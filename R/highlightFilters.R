#########################################################################/**
# @RdocFunction highlightFilters
#
# @title "Highlights data points in a plotted profile to evaluate filtering"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{object}{A QDNAseqReadCounts object ...}
#   \item{col}{...}
#   \item{blacklist}{...}
#   \item{mappability}{...}
#   \item{tgr}{...}
#   \item{bases}{...}
#   \item{...}{Further arguments to points.}
# }
#
# @author "IS"
#
# \seealso{
#   ...
# }
#
# @keyword IO
#*/#########################################################################
setMethod('highlightFilters', signature=c(object='QDNAseqReadCounts',
  col='character', blacklist='numeric', mappability='numeric', tgr='numeric',
  bases='numeric'), definition=function(object, col='red', blacklist=0,
  mappability=50, tgr=2, bases=100, ...) {
  condition <- condition <- rep(TRUE, times=nrow(object))
  if (!is.na(bases))
    condition <- condition & fData(object)$bases >= bases
  if (!is.na(blacklist))
    condition <- condition & fData(object)$blacklist <= blacklist
  if (!is.na(mappability))
    condition <- condition & fData(object)$mappability >= mappability
  if (!is.na(tgr))
    condition <- condition & !is.na(fData(object)$tgr) &
      abs(fData(object)$tgr) <= tgr*sd(fData(object)$tgr, na.rm=TRUE)
  all.chrom <- chromosomes(object)
  all.chrom.lengths <- aggregate(bpend(object),
    by=list(chromosome=all.chrom), max)
  chrom.lengths <- all.chrom.lengths$x
  names(chrom.lengths) <- all.chrom.lengths$chromosome
  chrom <- all.chrom[!condition]
  uni.chrom <- unique(chrom)
  chrom.lengths <- chrom.lengths[as.character(uni.chrom)]
  pos <- as.numeric(bpstart(object)[!condition])
  for (i in uni.chrom)
    pos[chrom > i] <- pos[chrom > i] + chrom.lengths[as.character(i)]
  copynumber <- copynumber(object)[!condition, , drop=FALSE]
  i <- 1
  pointcol <- col
  ylim <- par('usr')[3:4]
  points(pos, copynumber[, i], cex=.1, col=pointcol)
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
  chrom.in.plot <- unique(chromosomes(object)[binFilter(object)])
  has.data <- !is.na(copynumber(object)[, i])
  message('Highlighted ', format(sum(!condition[has.data &
    chromosomes(object) %in% chrom.in.plot]), big.mark=','), ' bins.')
})

# EOF

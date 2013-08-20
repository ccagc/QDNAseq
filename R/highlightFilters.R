#########################################################################/**
# @RdocFunction highlightFilters
#
# @alias highlightFilters,QDNAseqReadCounts-method
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
#   \item{mappability}{...}
#   \item{blacklist}{...}
#   \item{residual}{...}
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
setMethod('highlightFilters', signature=c(object='QDNAseqReadCounts'),
  definition=function(object, col='red', mappability=50, blacklist=0,
  residual=1, bases=100, ...) {
  condition <- condition <- rep(TRUE, times=nrow(object))
  if (!is.na(bases))
    condition <- condition & fData(object)$bases >= bases
  if (!is.na(blacklist))
    condition <- condition & fData(object)$blacklist <= blacklist
  if (!is.na(mappability))
    condition <- condition & fData(object)$mappability >= mappability
  if (!is.na(residual))
    condition <- condition & !is.na(fData(object)$residual) &
      abs(fData(object)$residual) <= residual*sd(fData(object)$residual,
      na.rm=TRUE)
  all.chrom <- chromosomes(object)
  all.chrom.lengths <- aggregate(bpend(object),
    by=list(chromosome=all.chrom), FUN=max)
  chrom.lengths <- all.chrom.lengths$x
  names(chrom.lengths) <- all.chrom.lengths$chromosome
  chrom <- all.chrom[!condition]
  uni.chrom <- unique(chrom)
  chrom.lengths <- chrom.lengths[as.character(uni.chrom)]
  pos <- as.numeric(bpstart(object)[!condition])
  for (i in uni.chrom)
    pos[chrom > i] <- pos[chrom > i] + chrom.lengths[as.character(i)]
  copynumber <- copynumber(object)[!condition, , drop=FALSE]
  if (ncol(object) > 1L)
    message('Multiple samples present in input, only using first sample: ',
      sampleNames(object)[1L])
  i <- 1
  pointcol <- col
  ylim <- par('usr')[3:4]
  points(pos, copynumber[, i], cex=0.1, col=pointcol)
  amps <- copynumber[, i]
  amps[amps < ylim[2]] <- NA_real_
  amps[!is.na(amps)] <- ylim[2] + 0.01 * (ylim[2]-ylim[1])
  dels <- copynumber[, i]
  dels[dels > ylim[1]] <- NA_real_
  dels[!is.na(dels)] <- ylim[1] - 0.01 * (ylim[2]-ylim[1])
  par(xpd=TRUE)
  points(pos, amps, pch=24, col=pointcol, bg=pointcol, cex=0.5)
  points(pos, dels, pch=25, col=pointcol, bg=pointcol, cex=0.5)
  par(xpd=FALSE)
  chrom.in.plot <- unique(chromosomes(object)[binFilter(object)])
  has.data <- !is.na(copynumber(object)[, i])
  num <- sum(!condition[has.data & chromosomes(object) %in% chrom.in.plot])
  message('Highlighted ', format(num, big.mark=','), ' bins.')
  return(invisible(num))
})

# EOF

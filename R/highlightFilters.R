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
#   \item{residual}{...}
#   \item{blacklist}{...}
#   \item{mappability}{...}
#   \item{bases}{...}
#   \item{type}{...}
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
  definition=function(object, col='red', residual=NA, blacklist=NA,
  mappability=NA, bases=NA, type='union', ...) {

  condition <- rep(TRUE, times=nrow(object))
  if (match.arg(type, c('union', 'intersection')) == 'intersection') {
    if (!is.na(residual)) {
      if (is.numeric(residual)) {
        condition <- condition & (is.na(fData(object)$residual) |
          abs(fData(object)$residual) > residual)
      } else if (residual) {
        condition <- condition & is.na(fData(object)$residual)
      }
    }
    if (!is.na(blacklist)) {
      if (is.numeric(blacklist)) {
        condition <- condition & fData(object)$blacklist > blacklist
      } else if (blacklist) {
        condition <- condition & fData(object)$blacklist != 0
      }
    }
    if (!is.na(mappability) && mappability != FALSE)
      condition <- condition & fData(object)$mappability < mappability
    if (!is.na(bases) && bases != FALSE)
      condition <- condition & fData(object)$bases < bases
  } else {
    if (!is.na(residual)) {
      if (is.numeric(residual)) {
        condition <- condition & !is.na(fData(object)$residual) &
          abs(fData(object)$residual) <= residual
      } else if (residual) {
        condition <- condition & !is.na(fData(object)$residual)
      }
    }
    if (!is.na(blacklist)) {
      if (is.numeric(blacklist)) {
        condition <- condition & fData(object)$blacklist <= blacklist
      } else if (blacklist) {
        condition <- condition & fData(object)$blacklist == 0
      }
    }
    if (!is.na(mappability) && mappability != FALSE)
      condition <- condition & fData(object)$mappability >= mappability
    if (!is.na(bases) && bases != FALSE)
      condition <- condition & fData(object)$bases >= bases
    condition <- !condition
  }

  i <- 1
  chrs <- unique(fData(object)$chromosome[binsToUse(object)])
  condition <- condition & fData(object)$bases > 0 &
    fData(object)$chromosome %in% chrs

  all.chrom <- chromosomes(object)
  all.chrom.lengths <- aggregate(bpend(object),
    by=list(chromosome=all.chrom), FUN=max)
  chrom.lengths <- all.chrom.lengths$x
  names(chrom.lengths) <- all.chrom.lengths$chromosome
  chrom <- all.chrom[condition]
  uni.chrom <- unique(chrom)
  chrom.lengths <- chrom.lengths[as.character(uni.chrom)]
  pos <- as.numeric(bpstart(object)[condition])
  for (j in uni.chrom)
    pos[chrom > j] <- pos[chrom > j] + chrom.lengths[as.character(j)]
  if ('copynumber' %in% assayDataElementNames(object)) {
    copynumber <- assayDataElement(object, 'copynumber')[condition, ,
      drop=FALSE]
  } else if ('corrected' %in% assayDataElementNames(object)) {
    copynumber <- assayDataElement(object, 'corrected')[condition, ,
      drop=FALSE]
  } else {
    copynumber <- assayDataElement(object, 'counts')[condition, , drop=FALSE]
  }
  if (ncol(object) > 1L)
    vmsg('Multiple samples present in input, only using first sample: ',
      sampleNames(object)[1L])
  pointcol <- col
  ylim <- par('usr')[3:4]
  points(pos, copynumber[, i], cex=0.1, col=pointcol)
  amps <- copynumber[, i]
  amps[amps <= ylim[2]] <- NA_real_
  amps[!is.na(amps)] <- ylim[2] + 0.01 * (ylim[2]-ylim[1])
  dels <- copynumber[, i]
  dels[dels >= ylim[1]] <- NA_real_
  dels[!is.na(dels)] <- ylim[1] - 0.01 * (ylim[2]-ylim[1])
  par(xpd=TRUE)
  points(pos, amps, pch=24, col=pointcol, bg=pointcol, cex=0.5)
  points(pos, dels, pch=25, col=pointcol, bg=pointcol, cex=0.5)
  par(xpd=FALSE)

  num <- sum(condition)
  vmsg('Highlighted ', format(num, big.mark=','), ' bins.')
  return(invisible(num))
})

# EOF

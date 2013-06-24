#########################################################################/**
# @RdocFunction correctBins
#
# @alias correctBins,QDNAseqReadCounts-method
#
# @title "Correct binned read counts for GC content and mappability"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{object}{...}
#   \item{span}{...}
#   \item{family}{...}
#   \item{adjustIncompletes}{...}
#   \item{keepCounts}{...}
#   \item{storeResiduals}{...}
#   \item{...}{Further aguments to loess.}
# }
#
# \value{
#   Returns a named @list with elements ...
# }
#
# @author "IS"
#
# \seealso{
#   Internally, @see "stats::loess" is used to fit the regression model.
# }
#*/#########################################################################

setMethod('correctBins', signature=c(object='QDNAseqReadCounts'),
  definition=function(object, span=0.65, family='symmetric',
  adjustIncompletes=TRUE, keepCounts=TRUE, storeResiduals=TRUE, ...) {
  counts <- assayDataElement(object, 'counts')
  if (adjustIncompletes) {
    counts <- counts / fData(object)$bases * 100
    counts[fData(object)$bases == 0] <- 0
  }
  message('Performing correction for GC content and mappability:')
  if (length(span) == 1L)
    span <- rep(span, times=ncol(counts))
  if (length(family) == 1L)
    family <- rep(family, times=ncol(counts))
  if (length(span) != ncol(counts))
    stop('Parameter span has to be either a single value or a vector the ',
      'same length as there are samples in object.')
  if (length(family) != ncol(counts))
    stop('Parameter family has to be either a single value or a vector the ',
      'same length as there are samples in object.')
  if ('filter' %in% colnames(fData(object))) {
    condition <- fData(object)$filter
  } else {
    condition <- rep(TRUE, times=nrow(object))
  }
  used.span <- rep(NA_real_, times=ncol(counts))
  used.family <- rep(NA_real_, times=ncol(counts))
  corrected <- matrix(nrow=nrow(counts), ncol=ncol(counts),
    dimnames=dimnames(counts))
  residuals <- matrix(nrow=nrow(counts), ncol=ncol(counts),
    dimnames=dimnames(counts))
  gc <- round(fData(object)$gc)
  mappability <- round(fData(object)$mappability)
  median.counts <- aggregate(counts[condition, ], by=list(gc=gc[condition],
    mappability=mappability[condition]), FUN=median)
  median.counts <- median.counts[!is.na(median.counts$gc), ]
  median.counts <- median.counts[!is.na(median.counts$mappability), ]
  rownames(median.counts) <- paste(median.counts$gc, '-',
    median.counts$mappability, sep='')
  # if (plotting) {
    # x <- min(median.counts$mappability):max(median.counts$mappability)
    # y <- min(median.counts$gc):max(median.counts$gc)
    # m <- matrix(nrow=length(x), ncol=length(y), dimnames=list(x, y))
    # d <- sub('-counts\\.(txt|RData)', '-correction', basename(countsfile))
    # if (!file.exists(d))
      # dir.create(d)
    # setwd(d)
  # }
  for (i in seq_len(ncol(counts))) {
    if (is.na(span[i]) && is.na(family[i])) {
      message('  Skipping correction for sample ', colnames(counts)[i], '...')
      next
    }
    message('  Using span=', span[i], '\tand family=', family[i],
      ',\tcorrecting sample ', colnames(counts)[i], '...')
    vals <- median.counts[, i+2L]
    corvals <- counts[, i]
    try({
      l <- loess(vals ~ median.counts$gc * median.counts$mappability,
        span=span[i], family=family[i]) # , ...)
      fit <- l$fitted
      names(fit) <- rownames(median.counts)
      residuals[, i] <- corvals - fit[paste(gc, '-', mappability, sep='')]
      correction <- median(fit, na.rm=TRUE) - fit
      corvals <- corvals + correction[paste(gc, '-', mappability, sep='')]
      corvals <- corvals - min(corvals, na.rm=TRUE)
      used.span[i] <- span[i]
      used.family[i] <- family[i]
      # if (plotting) {
        # pnga4(paste(colnames(counts)[i], '-counts.png', sep=''))
        # for (j in 1:nrow(median.counts))
          # m[as.character(median.counts[j, 'mappability']), as.character(median.counts[j, 'gc'])] <- median.counts[j, i+2]
        # image(x, y, m, col=paste('#', c(sprintf('%02X', 0:255), rep('FF', 256)), c(rep('FF', 256), sprintf('%02X', 255:0)), sprintf('%02X', 255), sep=''), xlab='mappability', ylab='GC content', main=paste('Median Read Counts -', colnames(counts)[i]), zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))))
        # contour(x, y, m, nlevels=20, zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))), add=TRUE)
        # dev.off()
        # pnga4(paste(colnames(counts)[i], '-fit.png', sep=''))
        # for (j in 1:nrow(median.counts))
          # m[as.character(median.counts[j, 'mappability']), as.character(median.counts[j, 'gc'])] <- fit[j]
        # image(x, y, m, col=paste('#', c(sprintf('%02X', 0:255), rep('FF', 256)), c(rep('FF', 256), sprintf('%02X', 255:0)), sprintf('%02X', 255), sep=''), xlab='mappability', ylab='GC content', main=paste('Loess Fit -', colnames(counts)[i]), zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))))
        # contour(x, y, m, nlevels=20, zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))), add=TRUE)
        # dev.off()
        # pnga4(paste(colnames(counts)[i], '-residuals.png', sep=''))
        # for (j in 1:nrow(median.counts))
          # m[as.character(median.counts[j, 'mappability']), as.character(median.counts[j, 'gc'])] <- l$residuals[j]
        # image(x, y, m, col=paste('#', c(sprintf('%02X', 0:255), rep('FF', 256)), c(rep('FF', 256), sprintf('%02X', 255:0)), sprintf('%02X', 255), sep=''), xlab='mappability', ylab='GC content', main=paste('Loess Residuals -', colnames(counts)[i]), zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))))
        # contour(x, y, m, nlevels=20, zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))), add=TRUE)
        # dev.off()
      # }
    }, silent=TRUE)
    corrected[, i] <- corvals
  }
  corrected <- round(corrected, digits=2L)
  # if (plotting) {
    # setwd('..')
    # pnga4(sub('-counts\\.(txt|RData)', '-residuals.png', basename(countsfile)))
    # image(x, y, m, col=paste('#', c(sprintf('%02X', 0:255), rep('FF', 256)), c(rep('FF', 256), sprintf('%02X', 255:0)), sprintf('%02X', 255), sep=''), xlab='mappability', ylab='GC content', main='Median Residuals', zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))))
    # contour(x, y, m, nlevels=20, zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))), add=TRUE)
    # dev.off()
  # }
  object$loess.span <- used.span
  object$loess.family <- used.family
  assayDataElement(object, 'corrected') <- corrected
  if (!keepCounts)
    assayDataElement(object, 'counts') <- NULL
  if (storeResiduals)
    assayDataElement(object, 'residuals') <- residuals
  message('Done.')
  object
})

# a subfunction for performing the correction on one sample at a time could
# be defined here.

# EOF

#########################################################################/**
# @RdocFunction segmentData
#
# @alias segmentData,QDNAseqReadCounts,logical-method
# @alias segmentData,cghRaw,missing-method
#
# @title "Segments and calls total copy numbers"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{object}{...}
#   \item{weights}{Either @TRUE or a vector of weights. If @TRUE,
#     loess residuals are used as weights.}
#   \item{...}{Additional arguments passed to @see "DNAcopy::segment".}
# }
#
# \value{
#   Returns ...
# }
#
# @author "IS"
#
# \seealso{
#   Internally, @see "DNAcopy::segment" of the \pkg{DNAcopy} package,
#   which implements the CBS method, is used to segment the data.
# }
#
#*/#########################################################################
setMethod('segmentData', signature=c(object='QDNAseqReadCounts'),
  definition=function(object, weights=TRUE, ...) {
  if (length(weights) == 1L & weights) {
    if ('residual' %in% colnames(fData(object))) {
      message('Using median loess residuals of control data set as ',
        'segmentation weights.')
      residual <- fData(object)$residual
    } else if ('residuals' %in% assayDataElementNames(object)) {
      message('Using median loess residuals as segmentation weights.')
      residual <- apply(assayDataElement(object, 'residuals'), 1, median,
        na.rm=TRUE)
    } else {
      stop('No loess residuals found. Please provide a vector of weights ',
        'or specify weights=FALSE.')
    }
    if (any(is.na(residual[binFilter(object)]))) {
      message('Filtering out ', sum(is.na(residual[binFilter(object)])),
        ' bins with missing residuals.')
      binFilter(object) <- binFilter(object) & !is.na(residual)
    }
    residual <- residual[binFilter(object)]
    residual <- abs(residual)
    residual[residual == 0] <- min(residual[residual != 0], na.rm=TRUE)
    weights <- 1/residual
  }
  copynumber <- copynumber(object)[binFilter(object), , drop=FALSE]
  CNA.object <- CNA(genomdat=copynumber, chrom=chromosomes(object)[binFilter(object)],
    maploc=bpstart(object)[binFilter(object)], data.type="logratio",
    sampleid=paste(sampleNames(object), ':', 1:ncol(object), 'of',
      ncol(object), sep=''))
  message('Start data segmentation ...')
  if (length(weights) == 1L && !weights) {
    segmented <- segment(CNA.object, ...)
  } else {
    segmented <- segment(CNA.object, weights=weights, ...)
  }
  numclone <- segmented$output$num.mark
  smrat <- segmented$output$seg
  numsmrat <- cbind(smrat, numclone)
  repdata <- function(row) {
    rep(row[1L], times=row[2L])
  }
  makelist <- apply(numsmrat, MARGIN=1L, FUN=repdata)
  joined <- unlist(makelist)
  ## Not needed anymore
  rm(list="makelist")
  joined <- matrix(joined, ncol=ncol(object), byrow=FALSE)
  rownames(joined) <- rownames(copynumber)
  colnames(joined) <- colnames(copynumber)
  segmented(object) <- joined
  object
})

setMethod('segmentData', signature=c(object='cghRaw'),
  definition=function(object, weights, ...) {
  CGHcall::segmentData(object, ...)
})

# EOF

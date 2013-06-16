#########################################################################/**
# @RdocFunction normalize
#
# @alias normalize,QDNAseqReadCounts-method
# @alias normalize,cghRaw-method
# 
# @title "Normalizes binned read counts"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{object}{...}
#   \item{method}{A @character string specifying ...}
#   \item{logTransform}{If @TRUE, ..., otherwise, ...}
#   \item{smoothOutliers}{If @TRUE, ..., otherwise, ...}
#   \item{...}{Further arguments to DNAcopy::smooth.CNA}
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
## Adapted from package CGHcall
setMethod('normalize', signature=c(object='QDNAseqReadCounts'),
  definition=function(object, method='median', smoothOutliers=TRUE,
  logTransform=TRUE, ...) {
  if ('filter' %in% colnames(fData(object))) {
    condition <- fData(object)$filter
  } else {
    condition <- rep(TRUE, times=nrow(object))
  }
  copynumber <- assayDataElement(object, 'corrected')[condition, , drop=FALSE]
  if (logTransform)
    copynumber <- log2(copynumber + 1)
  if (method == 'none') {
    message('Skipping normalization ...')
  } else {
    if (method == 'median') {
      message('Applying median normalization ...')
      ## TO DO: See matrixStats::rowMedians().
      values <- apply(copynumber, MARGIN=2L, FUN=median, na.rm=TRUE)
    } else if (method == 'mode') {
      message('Applying mode normalization ... ')
      values <- apply(copynumber, MARGIN=2L, FUN=function(x) { d <- density(x,
        na.rm=TRUE); d$x[which.max(d$y)] })
    }
    copynumber <- t(t(copynumber) - values)
  }
  if (smoothOutliers) {
    message('Smoothing outliers ...')
    CNA.object <- smooth.CNA(CNA(copynumber, fData(object)[condition,
      'chromosome'], fData(object)[condition, 'start'], data.type='logratio',
      presorted=TRUE), ...)
    copynumber <- as.matrix(CNA.object[, -(1:2), drop=FALSE])
  }
  copynumber2 <- matrix(nrow=nrow(object), ncol=ncol(object),
    dimnames=list(featureNames(object), sampleNames(object)))
  copynumber2[rownames(copynumber), ] <- copynumber
  assayDataElement(object, 'copynumber') <- copynumber2
  object
})

setMethod('normalize', signature=c(object='cghRaw'),
  definition=function(object, method='median', smoothOutliers=TRUE,
  logTransform=TRUE, ...) {
  CGHcall::normalize(input=object, method=method,
    smoothOutliers=smoothOutliers, ...)
})

# EOF

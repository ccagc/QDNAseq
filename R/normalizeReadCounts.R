#########################################################################/**
# @RdocFunction normalizeReadCounts
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
#   \item{obj}{...}
#   \item{method}{A @character string specifying ...}
#   \item{logTransform}{If @TRUE, ..., otherwise, ...}
#   \item{smoothOutliers}{If @TRUE, ..., otherwise, ...}
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
normalizeReadCounts <- function(obj, method='median', logTransform=TRUE,
  smoothOutliers=TRUE) {
  if ('filter' %in% colnames(fData(obj))) {
    condition <- fData(obj)$filter
  } else {
    condition <- rep(TRUE, times=nrow(obj))
  }
  copynumber <- assayDataElement(obj, 'corrected')[condition, , drop=FALSE]
  if (logTransform)
    copynumber <- log2(copynumber + 1)
  if (method == 'none') {
    message('Skipping normalization ...')
  } else {
    if (method == 'median') {
      message('Applying median normalization ...')
      # TO DO: See matrixStats::rowMedians().
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
    CNA.object <- smooth.CNA(CNA(copynumber, fData(obj)[condition,
      'chromosome'], fData(obj)[condition, 'start'], data.type='logratio',
      presorted=TRUE))
    copynumber <- as.matrix(CNA.object[,-(1:2), drop=FALSE])
  }
  copynumber2 <- matrix(nrow=nrow(obj), ncol=ncol(obj),
    dimnames=list(featureNames(obj), sampleNames(obj)))
  copynumber2[rownames(copynumber),] <- copynumber
  assayDataElement(obj, 'copynumber') <- copynumber2
  obj
}

# EOF

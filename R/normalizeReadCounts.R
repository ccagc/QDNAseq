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
normalizeReadCounts <- function(obj, method='median', logTransform=TRUE, smoothOutliers=TRUE) {
  if (exists('filter', obj)) {
    condition <- obj[['filter']]
  } else {
    condition <- rep(TRUE, times=nrow(obj[['bins']]))
  }
  bins <- obj[['bins']][condition,]
  copynumber <- obj[['corrected']][condition,]
  if (logTransform)
    copynumber <- log2(copynumber + 1)
  if (method == 'none') {
    cat('Skipping normalization ... \n')
  } else {
    if (method == 'median') {
      cat('Applying median normalization ... \n')
      # TO DO: See matrixStats::rowMedians().
      values <- apply(copynumber, MARGIN=2L, FUN=median, na.rm=TRUE)
    } else if (method == 'mode') {
      cat('Applying mode normalization ... \n')
      values <- apply(copynumber, MARGIN=2L, FUN=function(x) { d <- density(x, na.rm=TRUE); d$x[which.max(d$y)] })
    }
    copynumber <- t(t(copynumber) - values)
  }
  if (smoothOutliers) {
    cat('Smoothing outliers ... \n')
    CNA.object <- smooth.CNA(CNA(copynumber, bins$chromosome, bins$start, data.type='logratio', presorted=TRUE))
    copynumber <- as.matrix(CNA.object[,-(1:2)])
  }
  obj[['copynumber']] <- matrix(nrow=nrow(obj[['corrected']]), ncol=ncol(obj[['corrected']]), dimnames=dimnames(obj[['corrected']]))
  obj[['copynumber']][rownames(copynumber),] <- copynumber
  obj
}

# EOF

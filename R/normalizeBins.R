#########################################################################/**
# @RdocFunction normalizeBins
#
# @alias normalizeBins,QDNAseqReadCounts-method
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
## Adapted from CGHcall::normalize()
setMethod('normalizeBins', signature=c(object='QDNAseqReadCounts'),
  definition=function(object, method='median', smoothOutliers=TRUE,
  logTransform=TRUE, ...) {

  # Extract corrected counts
  counts <- assayDataElement(object, 'corrected')

  # Sanity check
  if (is.null(counts)) {
    stop(sprintf("Cannot normalize bins. %s object has no 'corrected' assay data.", class(object)[1L]))
  }

  # Extract annotation data
  fData <- fData(object);

  # Filter?
  if ('filter' %in% colnames(fData(object))) {
    keep <- fData$filter
    copynumber <- counts[keep, , drop=FALSE]
    fData <- fData[keep, , drop=FALSE];
  } else {
    copynumber <- counts
  }

  # Sanity check
  stopifnot(is.matrix(copynumber))

  # Log transform?
  if (logTransform)
    copynumber <- log2(copynumber + 1)

  if (method == 'none') {
    message('Skipping normalization ...')
  } else {
    if (method == 'median') {
      message('Applying median normalization ...')
      values <- colMedians(copynumber, na.rm=TRUE)
    } else if (method == 'mode') {
      message('Applying mode normalization ... ')
      values <- apply(copynumber, MARGIN=2L, FUN=function(x) {
        d <- density(x, na.rm=TRUE); d$x[which.max(d$y)]
      })
    }
    copynumber <- t(t(copynumber) - values)
  }

  # Smooth outliers?
  if (smoothOutliers) {
    message('Smoothing outliers ...')
    CNA.object <- CNA(copynumber, chrom=fData[,'chromosome'], maploc=fData[,'start'], data.type='logratio', presorted=TRUE)
    CNA.object <- smooth.CNA(CNA.object)
    CNA.object <- CNA.object[, -(1:2), drop=FALSE]
    copynumber <- as.matrix(CNA.object)
  }

  # Expand to full set of bins
  copynumber2 <- matrix(NA_real_, nrow=nrow(object), ncol=ncol(object),
               dimnames=list(featureNames(object), sampleNames(object)))
  copynumber2[rownames(copynumber), ] <- copynumber

  # Assign
  assayDataElement(object, 'copynumber') <- copynumber2

  object
})

# EOF

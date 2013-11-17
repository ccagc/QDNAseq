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
#   \item{force}{...}
#   \item{...}{Not used.}
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
  logTransform=TRUE, force=FALSE, ...) {

  if (!force && 'segmented' %in% assayDataElementNames(object))
    stop('Data has already been segmented. Changing the normalization will ',
      'remove segmentation (and possible calling) results. Please specify ',
      'force=TRUE, if you want this.')
  if ('segmented' %in% assayDataElementNames(object))
    assayDataElement(object, 'segmented') <- NULL
  if ('calls' %in% assayDataElementNames(object)) {
    assayDataElement(object, 'calls') <- NULL
    assayDataElement(object, 'probloss') <- NULL
    assayDataElement(object, 'probnorm') <- NULL
    assayDataElement(object, 'probgain') <- NULL
    if ('probdloss' %in% assayDataElementNames(object))
      assayDataElement(object, 'probdloss') <- NULL
    if ('probamp' %in% assayDataElementNames(object))
      assayDataElement(object, 'probamp') <- NULL
  }

  # Extract corrected counts
  copynumber <- assayDataElement(object, 'corrected')

  # Sanity check
  if (is.null(copynumber)) {
    stop(sprintf("Cannot normalize bins. %s object has no 'corrected' assay data.", class(object)[1L]))
  }

  # Extract annotation data
  fData <- fData(object)

  # Filter?
  condition <- binsToUse(object)

  # Sanity check
  stopifnot(is.matrix(copynumber))

  # Log transform?
  if (logTransform)
    copynumber <- log2(copynumber + 1)

  if (method == 'none') {
    vmsg('Skipping normalization ...')
  } else {
    if (method == 'median') {
      vmsg('Applying median normalization ...')
      values <- colMedians(copynumber[condition, , drop=FALSE], na.rm=TRUE)
    } else if (method == 'mode') {
      vmsg('Applying mode normalization ... ')
      values <- apply(copynumber[condition, , drop=FALSE], MARGIN=2L,
        FUN=function(x) {
        d <- density(x, na.rm=TRUE); d$x[which.max(d$y)]
      })
    }
    copynumber <- t(t(copynumber) - values)
  }

  # Smooth outliers?
  if (smoothOutliers) {
    vmsg('Smoothing outliers ...')
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

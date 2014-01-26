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
#   \item{object}{A @see "QDNAseqReadCounts" object with \code{corrected} data.}
#   \item{method}{A @character string specifying the normalization method.
#     Choices are "none", "mean", "median" (default), or "mode". A partial
#     string sufficient to uniquely identify the choice is permitted.}
#   \item{smoothOutliers}{If @TRUE (default), @see "DNAcopy::smooth.CNA" is
#     used to detect and smooth outliers.}
#   \item{logTransform}{If @TRUE (default), data will be log2-transformed.}
#   \item{logOffset}{An offset to be added prior to log2-transformation to
#     avoid non-positive numbers. Ignored if \code{logTransform} is @FALSE.}
#   \item{force}{Running this function will remove possible segmentation and
#     calling results. When they are present, running requires specifying
#     \code{force} is @TRUE.}
#   \item{...}{Further arguments to @see "DNAcopy::smooth.CNA".}
# }
#
# \value{
#   Returns a @see "QDNAseqReadCounts" object with the assay data element
#   \code{copynumbers} added.  If \code{logTransform} is @TRUE, these
#   signals are log2 transformed after adding the \code{logOffset} offset.
# }
#
# @author "IS"
#
#*/#########################################################################
## Adapted from CGHcall::normalize()
setMethod("normalizeBins", signature=c(object="QDNAseqReadCounts"),
  definition=function(object, method=c("median", "mean", "mode", "none"),
  smoothOutliers=TRUE, logTransform=TRUE, logOffset=2^-10, force=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'object':
  if (!force && ("segmented" %in% assayDataElementNames(object)))
    stop("Data has already been segmented. Changing the normalization will ",
      "remove segmentation (and possible calling) results. Please specify ",
      "force=TRUE, if you want this.")

  # Argument 'method':
  method <- match.arg(method);

  if ("segmented" %in% assayDataElementNames(object))
    assayDataElement(object, "segmented") <- NULL
  if ("calls" %in% assayDataElementNames(object)) {
    assayDataElement(object, "calls") <- NULL
    assayDataElement(object, "probloss") <- NULL
    assayDataElement(object, "probnorm") <- NULL
    assayDataElement(object, "probgain") <- NULL
    if ("probdloss" %in% assayDataElementNames(object))
      assayDataElement(object, "probdloss") <- NULL
    if ("probamp" %in% assayDataElementNames(object))
      assayDataElement(object, "probamp") <- NULL
  }

  # Extract corrected counts
  copynumber <- assayDataElement(object, "corrected")

  # Sanity check
  if (is.null(copynumber)) {
    stop("Cannot normalize bins. Please run correctBins() first.")
  }

  # Extract annotation data
  fData <- fData(object)

  # Sanity check
  stopifnot(is.matrix(copynumber))

  # Log transform?
  if (logTransform) {
    # remove negative values and add offset
    copynumber <- pmax(copynumber, 0)
    copynumber <- copynumber + logOffset
    copynumber <- log2(copynumber)
  }

  # Filter
  condition <- binsToUse(object)

  if (method == "none") {
    vmsg("Skipping normalization ...")
  } else {
    vmsg("Applying ", method, " normalization ...")
    if (method == "mean") {
      values <- colMeans(copynumber[condition, , drop=FALSE], na.rm=TRUE)
    } else if (method == "median") {
      values <- colMedians(copynumber[condition, , drop=FALSE], na.rm=TRUE)
    } else if (method == "mode") {
      values <- apply(copynumber[condition, , drop=FALSE], MARGIN=2L,
        FUN=function(x) {
        d <- density(x, na.rm=TRUE); d$x[which.max(d$y)]
      })
    }
    # Normalize by subtraction (log) or division (linear):
    if (logTransform) {
      copynumber <- scale(copynumber, center=values, scale=FALSE)
    } else {
      copynumber <- scale(copynumber, center=FALSE, scale=values)
    }
  }

  # Smooth outliers?
  if (smoothOutliers) {
    vmsg("Smoothing outliers ...")
    CNA.object <- CNA(copynumber, chrom=fData[,"chromosome"],
      maploc=fData[,"start"], data.type="logratio", presorted=TRUE)
    CNA.object <- smooth.CNA(CNA.object, ...)
    CNA.object <- CNA.object[, -(1:2), drop=FALSE]
    copynumber <- as.matrix(CNA.object)
    # Not needed anymore
    CNA.object <- NULL
  }

  # Expand to full set of bins
  copynumber2 <- matrix(NA_real_, nrow=nrow(object), ncol=ncol(object),
    dimnames=list(featureNames(object), sampleNames(object)))
  copynumber2[rownames(copynumber), ] <- copynumber

  # Not needed anymore
  copynumber <- NULL

  # Assign
  assayDataElement(object, "copynumber") <- copynumber2

  # Not needed anymore
  copynumber2 <- NULL

  object
})

# EOF

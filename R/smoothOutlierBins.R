#########################################################################/**
# @RdocFunction smoothOutlierBins
#
# @alias smoothOutlierBins,QDNAseqCopyNumbers-method
#
# @title "Smooth outlier bins after normalization"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{object}{A @see "QDNAseqCopyNumbers" object with \code{copynumber}
#         data.}
#     \item{logTransform}{If @TRUE (default), data will be log2-transformed.}
#     \item{force}{Running this function will remove possible segmentation and
#         calling results. When they are present, running requires specifying
#         \code{force} is @TRUE.}
#     \item{...}{Additional arguments passed to @see "DNAcopy::smooth.CNA".}
# }
#
# \value{
#     Returns a @see "QDNAseqCopyNumbers" object with the values for outliers
#     smoothed. See @see "DNAcopy::smooth.CNA" for more details. If
#     \code{logTransform} is @TRUE, these signals are log2-transformed prior
#     to smoothing, but afterwards back-transformed..
# }
#
# \examples{
# data(LGG150)
# readCounts <- LGG150
# readCountsFiltered <- applyFilters(readCounts)
# readCountsFiltered <- estimateCorrection(readCountsFiltered)
# copyNumbers <- correctBins(readCountsFiltered)
# copyNumbersNormalized <- normalizeBins(copyNumbers)
# copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
# }
#
# @author "IS"
#
# @keyword manip
#*/#########################################################################
## Adapted from CGHcall::normalize()
setMethod("smoothOutlierBins", signature=c(object="QDNAseqCopyNumbers"),
    definition=function(object, logTransform=TRUE, force=FALSE, ...) {
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument "object":
    if (!force && ("segmented" %in% assayDataElementNames(object)))
        stop("Data has already been segmented. Smoothing the outliers will ",
            "remove segmentation (and possible calling) results. ",
            "Please specify force=TRUE, if you want this.")

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
    copynumber <- assayDataElement(object, "copynumber")

    # Extract annotation data
    fData <- fData(object)

    # Sanity check
    stopifnot(is.matrix(copynumber))

    # Log transform?
    if (logTransform)
        copynumber <- log2adhoc(copynumber)

    # Filter
    condition <- binsToUse(object)

    vmsg("Smoothing outliers ...", appendLF=FALSE)
    CNA.object <- CNA(copynumber[condition, , drop=FALSE],
        chrom=fData[condition, "chromosome"],
        maploc=fData[condition, "start"],
        data.type="logratio", presorted=TRUE)
    CNA.object <- smooth.CNA(CNA.object, ...)
    CNA.object <- CNA.object[, -(1:2), drop=FALSE]
    copynumber <- as.matrix(CNA.object)
    # Not needed anymore
    CNA.object <- NULL

    # Log transform?
    if (logTransform)
        copynumber <- log2adhoc(copynumber, inv=TRUE)

    # Expand to full set of bins
    copynumber2 <- matrix(NA_real_, nrow=nrow(object), ncol=ncol(object),
        dimnames=list(featureNames(object), sampleNames(object)))
    copynumber2[condition, ] <- copynumber

    # Not needed anymore
    copynumber <- NULL

    # Assign
    assayDataElement(object, "copynumber") <- copynumber2

    # Not needed anymore
    copynumber2 <- NULL

    vmsg()
    object
})

# EOF

#########################################################################/**
# @RdocFunction normalizeBins
#
# @alias normalizeBins,QDNAseqCopyNumbers-method
#
# @title "Normalizes binned read counts"
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
#     \item{method}{A @character string specifying the normalization method.
#         Choices are "mean", "median" (default), or "mode". A partial
#         string sufficient to uniquely identify the choice is permitted.}
#     \item{force}{Running this function will remove possible segmentation and
#         calling results. When they are present, running requires specifying
#         \code{force} is @TRUE.}
#     \item{verbose}{If @TRUE, verbose messages are produced.}
# }
#
# \value{
#     Returns a @see "QDNAseqCopyNumbers" object with the assay data element
#     \code{copynumber} scaled with the chosen method.
# }
#
# \examples{
# data(LGG150)
# readCounts <- LGG150
# readCountsFiltered <- applyFilters(readCounts)
# readCountsFiltered <- estimateCorrection(readCountsFiltered)
# copyNumbers <- correctBins(readCountsFiltered)
# copyNumbersNormalized <- normalizeBins(copyNumbers)
# }
#
# @author "IS"
#
# @keyword manip
#*/#########################################################################
## Adapted from CGHcall::normalize()
setMethod("normalizeBins", signature=c(object="QDNAseqCopyNumbers"),
    definition=function(object, method=c("median", "mean", "mode"),
    force=FALSE,
    verbose=getOption("QDNAseq::verbose", TRUE)) {

    oopts <- options("QDNAseq::verbose"=verbose)
    on.exit(options(oopts))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument "object":
    if (!force && ("segmented" %in% assayDataElementNames(object)))
        stop("Data has already been segmented. Re-normalizing will ",
            "remove segmentation (and possible calling) results. ",
            "Please specify force=TRUE, if you want this.")

    # Argument "method":
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
    copynumber <- assayDataElement(object, "copynumber")

    # Extract annotation data
    fData <- fData(object)

    # Sanity check
    stopifnot(is.matrix(copynumber))

    # Filter
    condition <- binsToUse(object)

    vmsg("Applying ", method, " normalization ...", appendLF=FALSE)
    if (method == "mean") {
        values <- colMeans2(copynumber, rows=condition, na.rm=TRUE)
    } else if (method == "median") {
        values <- colMedians(copynumber, rows=condition, na.rm=TRUE)
    } else if (method == "mode") {
        values <- apply(copynumber[condition, , drop=FALSE], MARGIN=2L,
          FUN=function(x) {
            d <- density(x, na.rm=TRUE)
	    d$x[which.max(d$y)]
        })
    }
    vmsg()
    if (0 %in% values) {
        vmsg("These samples cannot be normalized (", method, "=0):\n",
            paste(sampleNames(object)[values == 0], collapse=", "))
        values[values == 0] <- 1
    }
    copynumber <- scale(copynumber, center=FALSE, scale=values)

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

#########################################################################/**
# @RdocFunction compareToReference
#
# @alias compareToReference,QDNAseqCopyNumbers,numeric-method
#
# @title "Divide binned read counts with those of reference samples"
#
# @synopsis
#
# \description{
#    @get "title".
# }
#
# \arguments{
#     \item{object}{An object of class @see "QDNAseqCopyNumbers".}
#     \item{references}{A numeric vector of indexes of the reference sample.
#         Must be the same length as there are samples in object. When @NA, the
#         sample will be kept as is. When @FALSE, the sample will be removed
#         from the output. As an example, object contains three samples: tumor1,
#         tumor2, and normal2. There is no reference for tumor1, but normal2 is
#         a matched normal sample from the same patient as tumor2. The keep
#         tumor1 as is, but to divide tumor2 with normal2, argument references
#         should be \code{c(NA, 3, FALSE)}.}
#     \item{force}{Whether to force the operation even when downstream data will
#         be lost.}
# }
#
# \value{
#     Returns a @see "QDNAseqCopyNumbers" object with the desired samples
#     divided by the signal of their reference samples.
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
# # Note: the following command will compare the sample to itself, which
# # does not really make sense:
# tumorVsNormal <- compareToReference(copyNumbersSmooth, c(1))
# }
#
# @author "IS"
#
# @keyword manip
#*/#########################################################################
setMethod("compareToReference", signature=c(object="QDNAseqCopyNumbers",
    references="numeric"),
    definition=function(object, references, force=FALSE) {

    if (!force && "segmented" %in% assayDataElementNames(object))
        stop("Data has already been segmented. Comparing to reference will ",
            "remove segmentation (and possible calling) ",
            "results. Please specify force=TRUE, if you want this.")
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

    if (length(references) != ncol(object))
        stop("Parameter references must be a vector of equal length as there ",
            "are samples in object.")
    for (i in seq_along(references)) {
        if (!is.na(references[i]) && references[i] != FALSE) {
            test <- assayDataElement(object, "copynumber")[, i]
            reference <- assayDataElement(object, "copynumber")[, references[i]]
            ratio <- test / reference
            ratio[reference == 0] <- 0
            assayDataElement(object, "copynumber")[, i] <- ratio
            sampleNames(object)[i] <- paste(sampleNames(object)[i], " vs. ",
                sampleNames(object)[references[i]], sep="")
            object$expected.variance[i] <- object$expected.variance[i] * 2
        }
    }
    toremove <- which(!references)
    if (length(toremove) > 0L)
        object <- object[, -toremove]
    object
})

# EOF

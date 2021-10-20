#########################################################################/**
# @RdocFunction normalizeSegmentedBins
#
# @alias normalizeSegmentedBins,QDNAseqCopyNumbers-method
#
# @title "Normalize segmented bins"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{object}{An object of class QDNAseqCopyNumbers.}
#     \item{inter}{The interval in which the function should search for the
#         normal level.}
#     \item{force}{Whether to force execution when it causes removal of
#         downstream calling results.}
# }
# \details{
#     This function recursively searches for the interval containing the
#     most segmented data, decreasing the interval length in each
#     recursion. The recursive search makes the post-segmentation
#     normalization robust against local maxima. This function is
#     particularly useful for profiles for which, after segmentation,
#     the 0-level does not coincide with many segments. It is more or
#     less harmless to other profiles. We advise to keep the search
#     interval (inter) small, in particular at the positive (gain) side
#     to avoid that the 0-level is set to a common gain level.
# }
# \value{
#     Returns an object of class QDNAseqCopyNumbers with re-normalized data.
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
# copyNumbersSegmented <- segmentBins(copyNumbersSmooth)
# copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
# }
#
# @author "IS"
#
# \seealso{
#     Internally, @see "CGHcall::postsegnormalize" of the \pkg{CGHcall} package
#     is used.
# }
#
# @keyword manip
#*/#########################################################################
setMethod("normalizeSegmentedBins", signature=c(object="QDNAseqCopyNumbers"),
    definition=function(object, inter=c(-0.1, 0.1), force=FALSE) {

    if (!force && "calls" %in% assayDataElementNames(object))
        stop("Data has already been called. Re-normalizing will ",
            "remove calling results. ",
            "Please specify force=TRUE, if you want this.")
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

    seg <- makeCgh(object, chromosomeReplacements="auto")
    postseg <- postsegnormalize(seg, inter=inter)
    copynumber(object) <- log2adhoc(copynumber(postseg), inv=TRUE)
    segmented(object) <- log2adhoc(segmented(postseg), inv=TRUE)
    object
})

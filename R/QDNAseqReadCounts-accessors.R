#########################################################################/**
# @RdocClass QDNAseqReadCounts
#
# @title "Container for QDNAseq read count data"
#
# \description{
#     @get "title"
# }
#
# \section{Assay data elements}{
#     An object of this class contains (a subset) the following elements:
#     \describe{
#         \item{\code{counts}}{(@numeric) Binned read counts as non-negative
#             integers in \eqn{\{0,1,2,...\}}. An object with only this field is
#             created by @see "binReadCounts".}
#         \item{\code{fit}}{(@numeric; optional) Loess fit of "count" signals as
#             doubles. Normally, these should all be positive values, but a
#             small number of edge case bins might contain negatives, especially
#             when fitting unfiltered data. This element is added after calling
#             @see "estimateCorrection".}
#     }
# }
#
# \section{Missing values}{
#     The bin data (assay data) may contain missing values.
# }
#
#
# @author "IS"
#*/#########################################################################

setMethod('initialize', 'QDNAseqReadCounts',
    function(.Object, bins, counts, phenodata, ...) {
    
    if (inherits(bins, "data.frame"))
        bins <- AnnotatedDataFrame(bins)
    if (inherits(phenodata, "data.frame"))
        phenodata <- AnnotatedDataFrame(phenodata)
    callNextMethod(.Object, featureData=bins,
        assayData=assayDataNew(counts=counts),
        phenoData=phenodata, ...)
})

setMethod('counts', signature=c(object='QDNAseqReadCounts'),
    definition=function(object) {
    
    assayDataElement(object, 'counts')
})

setMethod('fit', signature=c(object='QDNAseqReadCounts'),
    definition=function(object) {
    
    assayDataElement(object, 'fit')
})

# setReplaceMethod('counts', signature=c(object='QDNAseqReadCounts',
#     value='matrix'), definition=function(object, value) {
#     if (nrow(value) == nrow(object)) {
#         assayDataElementReplace(object, 'counts', value)
#     } else {
#         value2 <- matrix(nrow=nrow(object), ncol=ncol(object),
#             dimnames=list(featureNames(object), sampleNames(object)))
#         value2[rownames(value), ] <- value
#         assayDataElementReplace(object, 'counts', value2)
#     }
# })

# setReplaceMethod('fit', signature=c(object='QDNAseqReadCounts',
#     value='matrix'), definition=function(object, value) {
#     if (nrow(value) == nrow(object)) {
#         assayDataElementReplace(object, 'fit', value)
#     } else {
#         value2 <- matrix(nrow=nrow(object), ncol=ncol(object),
#             dimnames=list(featureNames(object), sampleNames(object)))
#         value2[rownames(value), ] <- value
#         assayDataElementReplace(object, 'fit', value2)
#     }
# })

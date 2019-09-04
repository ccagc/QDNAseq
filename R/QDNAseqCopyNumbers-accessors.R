#########################################################################/**
# @RdocClass QDNAseqCopyNumbers
#
# @title "Container for QDNAseq read count data"
#
# \description{
#     @get "title"
# }
#
# \section{Assay data elements}{
#     An object of this class contains the following elements:
#     \describe{
#         \item{\code{copynumber}}{(@numeric) Corrected "count" signals in
#             \eqn{[0,+\infty)} An object with only this field is created by
#             @see "correctBins".}
#         \item{\code{segmented}}{(@numeric; optional) Segmented data in
#             \eqn{[0,+\infty)}, added by calling @see "segmentBins".}
#         \item{\code{calls}}{(@integer; optional) Calls as -2=deletion,
#             -1=loss, 0=normal, 1=gain, 2=amplification, added by calling
#             @see "callBins".}
#         \item{\code{probdloss}}{(@numeric; optional) Probabilities of
#             deletions in \eqn{[0,1]}, added by calling @see "callBins".}
#         \item{\code{probloss}}{(@numeric; optional) Probabilities of losses in
#             \eqn{[0,1]}, added by calling @see "callBins".}
#         \item{\code{probnorm}}{(@numeric; optional) Probabilities of normal
#             copy number in \eqn{[0,1]}, added by calling @see "callBins".}
#         \item{\code{probgain}}{(@numeric; optional) Probabilities of gains in
#             \eqn{[0,1]}, added by calling @see "callBins".}
#         \item{\code{probamp}}{(@numeric; optional) Probabilities of
#             amplifications in \eqn{[0,1]}, added by calling @see "callBins".}
#     }
# }
#
# \section{Missing values}{
#     The bin data (assay data) may contain missing values.
# }
#
# @author "IS"
#*/#########################################################################

setMethod('initialize', 'QDNAseqCopyNumbers',
    function(.Object, bins, copynumber, phenodata, ...) {
    
    if (class(bins) == 'data.frame')
        bins <- AnnotatedDataFrame(bins)
    if (class(phenodata) == 'data.frame')
        phenodata <- AnnotatedDataFrame(phenodata)
    callNextMethod(.Object, featureData=bins,
        assayData=assayDataNew(copynumber=copynumber),
        phenoData=phenodata, ...)
})

setMethod('copynumber', signature=c(object='QDNAseqCopyNumbers'),
    definition=function(object) {
    
    assayDataElement(object, 'copynumber')
})

setMethod('segmented', signature=c(object='QDNAseqCopyNumbers'),
    definition=function(object) {
    
    assayDataElement(object, 'segmented')
})

setMethod('calls', signature=c(object='QDNAseqCopyNumbers'),
    definition=function(object) {
    
    assayDataElement(object, 'calls')
})

setMethod('probdloss', signature=c(object='QDNAseqCopyNumbers'),
    definition=function(object) {
    
    assayDataElement(object, 'probdloss')
})

setMethod('probloss', signature=c(object='QDNAseqCopyNumbers'),
    definition=function(object) {
    
    assayDataElement(object, 'probloss')
})

setMethod('probnorm', signature=c(object='QDNAseqCopyNumbers'),
    definition=function(object) {
    
    assayDataElement(object, 'probnorm')
})

setMethod('probgain', signature=c(object='QDNAseqCopyNumbers'),
    definition=function(object) {
    
    assayDataElement(object, 'probgain')
})

setMethod('probamp', signature=c(object='QDNAseqCopyNumbers'),
    definition=function(object) {
    
    assayDataElement(object, 'probamp')
})

setReplaceMethod('copynumber',
    signature=c(object='QDNAseqCopyNumbers', value='matrix'),
    definition=function(object, value) {
    
    if (nrow(value) == nrow(object)) {
        assayDataElementReplace(object, 'copynumber', value)
    } else {
        value2 <- matrix(NA_real_, nrow=nrow(object), ncol=ncol(object),
            dimnames=list(featureNames(object), sampleNames(object)))
        value2[rownames(value), ] <- value
        assayDataElementReplace(object, 'copynumber', value2)
    }
})

setReplaceMethod('segmented',
    signature=c(object='QDNAseqCopyNumbers', value='matrix'),
    definition=function(object, value) {
    
    if (nrow(value) == nrow(object)) {
        assayDataElementReplace(object, 'segmented', value)
    } else {
        value2 <- matrix(NA_real_, nrow=nrow(object), ncol=ncol(object),
            dimnames=list(featureNames(object), sampleNames(object)))
        value2[rownames(value), ] <- value
        assayDataElementReplace(object, 'segmented', value2)
    }
})

setReplaceMethod('calls',
    signature=c(object='QDNAseqCopyNumbers', value='matrix'),
    definition=function(object, value) {
    
    if (nrow(value) == nrow(object)) {
        assayDataElementReplace(object, 'calls', value)
    } else {
        value2 <- matrix(NA_real_, nrow=nrow(object), ncol=ncol(object),
            dimnames=list(featureNames(object), sampleNames(object)))
        value2[rownames(value), ] <- value
        assayDataElementReplace(object, 'calls', value2)
    }
})

setReplaceMethod('probdloss',
    signature=c(object='QDNAseqCopyNumbers', value='matrix'),
    definition=function(object, value) {
    
    if (nrow(value) == nrow(object)) {
        assayDataElementReplace(object, 'probdloss', value)
    } else {
        value2 <- matrix(NA_real_, nrow=nrow(object), ncol=ncol(object),
            dimnames=list(featureNames(object), sampleNames(object)))
        value2[rownames(value), ] <- value
        assayDataElementReplace(object, 'probdloss', value2)
    }
})

setReplaceMethod('probloss',
    signature=c(object='QDNAseqCopyNumbers', value='matrix'),
    definition=function(object, value) {
    
    if (nrow(value) == nrow(object)) {
        assayDataElementReplace(object, 'probloss', value)
    } else {
        value2 <- matrix(NA_real_, nrow=nrow(object), ncol=ncol(object),
            dimnames=list(featureNames(object), sampleNames(object)))
        value2[rownames(value), ] <- value
        assayDataElementReplace(object, 'probloss', value2)
    }
})

setReplaceMethod('probnorm',
    signature=c(object='QDNAseqCopyNumbers', value='matrix'),
    definition=function(object, value) {
    
    if (nrow(value) == nrow(object)) {
        assayDataElementReplace(object, 'probnorm', value)
    } else {
        value2 <- matrix(NA_real_, nrow=nrow(object), ncol=ncol(object),
            dimnames=list(featureNames(object), sampleNames(object)))
        value2[rownames(value), ] <- value
        assayDataElementReplace(object, 'probnorm', value2)
    }
})

setReplaceMethod('probgain',
    signature=c(object='QDNAseqCopyNumbers', value='matrix'),
    definition=function(object, value) {
    
    if (nrow(value) == nrow(object)) {
        assayDataElementReplace(object, 'probgain', value)
    } else {
        value2 <- matrix(NA_real_, nrow=nrow(object), ncol=ncol(object),
            dimnames=list(featureNames(object), sampleNames(object)))
        value2[rownames(value), ] <- value
        assayDataElementReplace(object, 'probgain', value2)
    }
})

setReplaceMethod('probamp',
    signature=c(object='QDNAseqCopyNumbers', value='matrix'),
    definition=function(object, value) {
    
    if (nrow(value) == nrow(object)) {
        assayDataElementReplace(object, 'probamp', value)
    } else {
        value2 <- matrix(NA_real_, nrow=nrow(object), ncol=ncol(object),
            dimnames=list(featureNames(object), sampleNames(object)))
        value2[rownames(value), ] <- value
        assayDataElementReplace(object, 'probamp', value2)
    }
})

# EOF

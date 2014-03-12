#########################################################################/**
# @RdocFunction makeCgh
#
# @alias makeCgh,QDNAseqCopyNumbers-method
#
# @title "Constructs a 'cghRaw', 'cghSeg', or 'cghCall' object"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{object}{A @see "QDNAseqCopyNumbers" object.}
#     \item{filter}{If @TRUE, bins are filtered, otherwise not.}
#     \item{...}{Not used.}
# }
#
# \value{
#     Returns a @see "CGHbase::cghRaw" if the object has not been segmented,
#     a @see "CGHbase::cghSeg" if it has been segmented but not called,
#     or @see "CGHbase::cghCall" if it has been called as well.
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
# cgh <- makeCgh(copyNumbersSmooth)
# }
#
# @author "IS"
#
# @keyword manip
#*/#########################################################################
setMethod('makeCgh', signature=c(object='QDNAseqCopyNumbers'),
    definition=function(object, filter=TRUE, ...) {

    # Decide which cgh* class to create
    names <- assayDataElementNames(object)
    if ('calls' %in% names) {
        className <- 'cghCall'
        segmented(object) <- log2adhoc(segmented(object))
        copynumber(object) <- log2adhoc(copynumber(object))
    } else if ('segmented' %in% names) {
        className <- 'cghSeg'
        segmented(object) <- log2adhoc(segmented(object))
        copynumber(object) <- log2adhoc(copynumber(object))
    } else if ('copynumber' %in% names) {
        className <- 'cghRaw'
        copynumber(object) <- log2adhoc(copynumber(object))
    } else {
        stop("Cannot create a CGHbase::cgh* object without assay data element",
            " 'copynumber': ", paste(sQuote(names), collapse=", "))
    }

    # Filter bins?
    if (filter) {
        keep <- binsToUse(object)
        object <- object[keep,]
        keep <- NULL # Not needed anymore
    }

    # Coerce chromosomes to integer indices
    fData(object)$chromosome <- chromosomes(object)

    # Update column names
    names <- colnames(fData(object))
    names[names == 'chromosome'] <- 'Chromosome'
    names[names == 'start'] <- 'Start'
    names[names == 'end'] <- 'End'
    colnames(fData(object)) <- names

    # Instantiate choose cgh* object
    new(className, assayData=assayData(object),
        featureData=featureData(object),
        phenoData=phenoData(object))
})


# Methods for coercing a QDNAseqCopyNumbers object into
# cghRaw, cghSeg and cghCall object using as(from, to).
setAs("QDNAseqCopyNumbers", "cghRaw", function(from) {
    makeCgh(from, filter=FALSE)
})

setAs("QDNAseqCopyNumbers", "cghSeg", function(from) {
    makeCgh(from, filter=FALSE)
})

setAs("QDNAseqCopyNumbers", "cghCall", function(from) {
    makeCgh(from, filter=FALSE)
})

# EOF

#########################################################################/**
# @RdocFunction segmentBins
#
# @alias segmentBins,QDNAseqCopyNumbers-method
#
# @title "Segments normalized copy number data"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{object}{An object of class QDNAseqCopyNumbers.}
#     \item{smoothBy}{An optional integer value to perform smoothing before
#         segmentation by taking the mean of every smoothBy bins, and then
#         segment those means. Default is to perform no smoothing.}
#     \item{alpha}{Significance levels for the test to accept change-points.
#         Default is 1e-10.}
#     \item{undo.splits}{A character string specifying how change-points are to
#         be undone, if at all. Default is "sdundo", which undoes splits that
#         are not at least this many SDs apart. Other choices are
#         "prune", which uses a sum of squares criterion, and "none".}
#     \item{undo.SD}{The number of SDs between means to keep a split if
#         undo.splits="sdundo". Default is 1.0.}
#     \item{force}{Whether to force execution when it causes removal of
#         downstream calling results.}
#     \item{transformFun}{A function to transform the data with. This can be
#         the default "log2" for log2(x + .Machine$double.xmin),
#         "sqrt" for the Anscombe transform of sqrt(x * 3/8) which
#         stabilizes the variance, "none" for no transformation, or any
#         R function that performs the desired transformation and also its
#         inverse when called with parameter \code{inv=TRUE}.}
#     \item{...}{Additional arguments passed to @see "DNAcopy::segment".}
# }
#
# \value{
#     Returns an object of class QDNAseqCopyNumbers with segmentation results
#         added.
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
# }
#
# @author "IS"
#
# \seealso{
#     Internally, @see "DNAcopy::segment" of the \pkg{DNAcopy} package,
#     which implements the CBS method, is used to segment the data.
# }
#
# @keyword manip
# @keyword smooth
#*/#########################################################################
setMethod("segmentBins", signature=c(object="QDNAseqCopyNumbers"),
    definition=function(object, smoothBy=FALSE, alpha=1e-10,
    undo.splits="sdundo", undo.SD=1.0, force=FALSE,
    transformFun="log2", ... ) {

    if (!force && "calls" %in% assayDataElementNames(object))
        stop("Data has already been called. Re-segmentation will ",
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
    condition <- binsToUse(object)

    if (!is.numeric(smoothBy) || smoothBy <= 1) {
        vmsg("Performing segmentation:")
    } else {
        vmsg("Performing segmentation with smoothing over ",
            smoothBy, " bins:")
    }

    copynumber <- copynumber(object)
    copynumber[!condition, ] <- NA_real_
    if (is.character(transformFun)) {
        transformFun <- match.arg(transformFun, c("log2", "sqrt", "none"))
        if (transformFun == "log2") {
            transformFun <- QDNAseq:::log2adhoc
        } else if (transformFun == "sqrt") {
            transformFun <- QDNAseq:::sqrtadhoc
        }
    }
    if (!is.character(transformFun) || transformFun != "none") {
        transformFun <- match.fun(transformFun)
        copynumber <- transformFun(copynumber)
    }
    segmented <- matrix(NA_real_, nrow=nrow(copynumber), ncol=ncol(copynumber),
        dimnames=dimnames(copynumber))

    ## loop through samples
    for (s in seq_len(ncol(copynumber))) {
        vmsg("    Segmenting: ", sampleNames(object)[s],
            " (", s, " of ", ncol(object), ") ...", appendLF=FALSE)

        ## loop through chromosomes
        for (chr in unique(fData(object)$chromosome[condition])) {
            index <- fData(object)$chromosome == chr
            chrStarts <- fData(object)$start[index]
            ## smooth if needed
            if (is.na(smoothBy) || !smoothBy || smoothBy <= 1) {
                chrCopynumber <- copynumber[index, s]
            } else {
                binToBin <- 0:(sum(index)-1) %/% smoothBy
                chrCopynumber <- aggregate(copynumber[index, s],
                    by=list(binToBin), mean, na.rm=TRUE)$x
                chrStarts <- chrStarts[seq(from=1, by=smoothBy,
                    length.out=length(chrCopynumber))]
            }

            ## segment
            cna <- CNA(genomdat=chrCopynumber,
                chrom=chr, maploc=chrStarts, data.type="logratio",
                presorted=TRUE)
            segments <- segment(cna, verbose=0,
                alpha=alpha, undo.splits=undo.splits, undo.SD=undo.SD, ...)
            chrSegmented <- rep(NA_real_, length=length(chrCopynumber))
            for (i in 1:nrow(segments$output))
                chrSegmented[segments$segRows$startRow[i]:
                    segments$segRows$endRow[i]] <- segments$output$seg.mean[i]

            ## process results whether smoothed or not
            if (is.na(smoothBy) || !smoothBy || smoothBy <= 1) {
                segmented[index, s] <- chrSegmented
            } else {
                segmented[index, s] <- rep(chrSegmented, times=table(binToBin))
            }
        }
        vmsg()
    }
    segmented[is.na(copynumber)] <- NA_real_
    if (!is.character(transformFun) || transformFun != "none")
        segmented <- transformFun(segmented, inv=TRUE)
    segmented(object) <- segmented
    object
})

# EOF

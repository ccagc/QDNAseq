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
#         segment those means. Default (@FALSE) is to perform no smoothing.
#         \code{smoothBy=1L} is a special case that will not perform smoothing,
#         but will split the segmentation process by chromosome instead of by
#         sample. This has an effect when using parallel computing or specifying
#         seeds for random number generation.}
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
#     \item{seeds}{An optional integer vector of seeds for random number
#         generation, recycled as needed. Normally, the segmentation process is
#         split by sample, and provided seeds also used per sample. But when
#         smoothing is performed (or in the the special case of
#         \code{smoothBy=1L}), the process is split by chromosome, seeds used
#         per chromosome, and results not necessarily reproducible across
#         samples.}
#%     \item{segmentStatistic}{A character vector specifying which segment 
#%         statistic to use.}
#%     \item{storeSegmentObjects}{A boolean to indicate whether to store the raw
#%         DNAcopy objects within the QDNAseq objects. Segment objects can be
#%         retrieved as segmentObject in assayData 
#%         eg. "assayDataElement(object, 'segmentObject')"}
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
    transformFun="log2", seeds=NULL, 
    segmentStatistic="seg.mean", storeSegmentObjects=FALSE,
    ...) {

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

    copynumber <- assayDataElement(object, "copynumber")
    copynumber[!condition, ] <- NA_real_
    if (is.character(transformFun)) {
        transformFun <- match.arg(transformFun, c("log2", "sqrt", "none"))
        if (transformFun == "log2") {
            transformFun <- log2adhoc
        } else if (transformFun == "sqrt") {
            transformFun <- sqrtadhoc
        }
    }
    if (!is.character(transformFun) || transformFun != "none") {
        transformFun <- match.fun(transformFun)
        copynumber <- transformFun(copynumber)
    }

    # if (is.na(smoothBy) || !smoothBy || smoothBy <= 1) {
    if (is.na(smoothBy) || !smoothBy) {
        ## create a list of CNA objects that can be analyzed with *lapply()
        cna <- lapply(sampleNames(object), function(x)
            CNA(genomdat=copynumber[condition, x, drop=FALSE],
                chrom=factor(fData(object)$chromosome[condition],
                    levels=unique(fData(object)$chromosome), ordered=TRUE),
                maploc=fData(object)$start[condition], data.type="logratio",
                sampleid=x, presorted=TRUE))
        ## create a vector of messages to be printed
        msgs <- paste0("    Segmenting: ", sampleNames(object),
            " (", 1:ncol(object), " of ", ncol(object), ") ...")
        ## use sample names for indexing, they are available in the CNA objects
        ## as the name of the third column
        names(msgs) <- sampleNames(object)
        segments <- flapply(cna, FUN=function(x, ..., seeds=NULL) {
            vmsg(msgs[colnames(x)[3]])
            segment(x, alpha=alpha, undo.splits=undo.splits,
                    undo.SD=undo.SD, verbose=0, ...)
        }, ..., seeds=seeds)
        
        if(storeSegmentObjects)
          assayDataElementReplace(object, "segmentObj", segments) -> object
        
        segmentStatisticCol <- grep(segmentStatistic, 
                                    colnames(segments.summary(segments[[1]])))
        
        segmented <- do.call(cbind, lapply(segments, function(x)
            rep(segments.summary(x)[,segmentStatisticCol], x$output$num.mark)))
        dimnames(segmented) <- dimnames(copynumber[condition, , drop=FALSE])
    } else {
        cna <- lapply(unique(fData(object)$chromosome[condition]),
            function(chr) {
                index <- fData(object)$chromosome == chr
                chrStarts <- fData(object)$start[index]
                ## smooth if needed
                if (is.na(smoothBy) || !smoothBy || smoothBy <= 1) {
                    chrCopynumber <- copynumber[index, , drop=FALSE]
                } else {
                    binToBin <- 0:(sum(index)-1) %/% smoothBy
                    chrCopynumber <- aggregate(copynumber[index, , drop=FALSE],
                        by=list(binToBin=binToBin), mean, na.rm=TRUE)
                    chrCopynumber$binToBin <- NULL
                    chrCopynumber <- as.matrix(chrCopynumber)
                    chrStarts <- chrStarts[seq(from=1, by=smoothBy,
                        length.out=nrow(chrCopynumber))]
                }
                CNA(genomdat=chrCopynumber, chrom=chr, maploc=chrStarts,
                    data.type="logratio", sampleid=sampleNames(object),
                    presorted=TRUE)
            }
        )

        segments <- flapply(cna, FUN=function(x, ..., seeds=NULL) {
            vmsg("    Segmenting chromosome ", x$chrom[1], " ...")
            segment(x, alpha=alpha, undo.splits=undo.splits,
                    undo.SD=undo.SD, verbose=0, ...)
        }, ..., seeds=seeds)
        
        if(storeSegmentObjects)
          assayDataElementReplace(object, "segmentObj", segments) -> object
        
        segmentStatisticCol <- grep(segmentStatistic, 
                                    colnames(segments.summary(segments[[1]])))
        
        segmented <- do.call(rbind, lapply(segments, function(x) {
            chrSegmented <- matrix(NA_real_, nrow=nrow(x$data),
                ncol=ncol(x$data)-2,
                dimnames=list(NULL, colnames(x$data)[-(1:2)]))
            for (i in 1:nrow(x$output))
                chrSegmented[x$segRows$startRow[i]:x$segRows$endRow[i],
                    x$output$ID[i]] <- segments.summary(x)[i, segmentStatisticCol]
            ## process results whether smoothed or not
            index <- fData(object)$chromosome == x$data$chrom[1]
            binToBin <- 0:(sum(index)-1) %/% smoothBy
            if (is.na(smoothBy) || !smoothBy || smoothBy <= 1) {
                # chrSegmented
            } else {
                chrSegmented <- apply(chrSegmented, 2, rep,
                    times=table(binToBin))
            }
            rownames(chrSegmented) <- rownames(copynumber)[index]
            chrSegmented
        }))
    }
    if (!is.character(transformFun) || transformFun != "none")
        segmented <- transformFun(segmented, inv=TRUE)
    segmented(object) <- segmented
    segmented(object)[is.na(copynumber)] <- NA_real_
    object
})

# EOF

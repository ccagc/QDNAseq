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
#         sample.}
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
#%     \item{segmentStatistic}{A character vector specifying which segment 
#%         statistic to use.}
#%     \item{storeSegmentObjects}{A boolean to indicate whether to store the raw
#%         DNAcopy objects within the QDNAseq objects. Segment objects can be
#%         retrieved as segmentObject in assayData 
#%         eg. "assayDataElement(object, 'segmentObject')"}
#     \item{...}{Additional arguments passed to @see "DNAcopy::segment".}
#%     \item{verbose}{If @TRUE, verbose messages are produced.}
# }
#
# \value{
#     Returns an object of class QDNAseqCopyNumbers with segmentation results
#         added.
# }
#
# \section{Numerical reproducibility}{
#  This method make use of random number generation (RNG) via the
#  @see "DNAcopy::segment" used internally.  Because of this, calling the
#  method with the same input data multiple times will each time give slightly
#  different results.  To get numerically reproducible results, the random
#  seed must be fixed, e.g. by using `set.seed()` at the top of the script.
# }
#
# \section{Parallel processing}{
#   This function uses \pkg{future} to segment samples in parallel.
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
#     which implements the CBS method [1,2], is used to segment the data.
# }
#
# \references{
#  [1] A.B. Olshen, E.S. Venkatraman (aka Venkatraman E. Seshan), R. Lucito
#      and M. Wigler, \emph{Circular binary segmentation for the analysis of
#      array-based DNA copy number data}, Biostatistics, 2004 \cr
#  [2] E.S. Venkatraman and A.B. Olshen, \emph{A faster circular binary
#      segmentation algorithm for the analysis of array CGH data},
#      Bioinformatics, 2007 \cr
# }
#
# @keyword manip
# @keyword smooth
#*/#########################################################################
setMethod("segmentBins", signature=c(object="QDNAseqCopyNumbers"),
    definition=function(object, smoothBy=FALSE, alpha=1e-10,
    undo.splits="sdundo", undo.SD=1.0, force=FALSE,
    transformFun="log2",
    segmentStatistic="seg.mean", storeSegmentObjects=FALSE,
    ..., verbose=getOption("QDNAseq::verbose", TRUE)) {

    if ("seeds" %in% names(list(...))) {
      .Defunct("Argument 'seeds' (integer) is no longer supported and ignored.")
    }
    
    oopts <- options("QDNAseq::verbose"=verbose)
    on.exit(options(oopts))

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
        cna <- lapply(sampleNames(object), FUN=function(x) {
            CNA(genomdat=copynumber[condition, x, drop=FALSE],
                chrom=factor(fData(object)$chromosome[condition],
                    levels=unique(fData(object)$chromosome), ordered=TRUE),
                maploc=fData(object)$start[condition], data.type="logratio",
                sampleid=x, presorted=TRUE)
        })
        ## create a vector of messages to be printed
        msgs <- paste0("    Segmenting: ", sampleNames(object),
            " (", 1:ncol(object), " of ", ncol(object), ") ...")
        segments <- future_mapply(cna, msgs, FUN=function(x, msg, ...) {
            vmsg(msg)
            segment(x, alpha=alpha, undo.splits=undo.splits,
                    undo.SD=undo.SD, verbose=0, ...)
        }, MoreArgs = list(...), SIMPLIFY = FALSE, USE.NAMES = FALSE, future.seed=TRUE)
        
        if(storeSegmentObjects)
          object <- assayDataElementReplace(object, "segmentObj", segments)
        
        segmentStatisticCol <- grep(segmentStatistic, 
                                    colnames(segments.summary(segments[[1]])))
        
        segmented <- do.call(cbind, lapply(segments, FUN=function(x)
            rep(segments.summary(x)[,segmentStatisticCol], times=x$output$num.mark)))
        dimnames(segmented) <- dimnames(copynumber[condition, , drop=FALSE])
    } else {
        cna <- lapply(unique(fData(object)$chromosome[condition]),
            FUN=function(chr) {
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

        segments <- future_lapply(cna, FUN=function(x, ...) {
            vmsg("    Segmenting chromosome ", x$chrom[1], " ...")
            segment(x, alpha=alpha, undo.splits=undo.splits,
                    undo.SD=undo.SD, verbose=0, ...)
        }, ..., future.seed=TRUE)
        
        if(storeSegmentObjects)
          object <- assayDataElementReplace(object, "segmentObj", segments)
        
        segmentStatisticCol <- grep(segmentStatistic, 
                                    colnames(segments.summary(segments[[1]])))
        
        segmented <- do.call(rbind, args = lapply(segments, FUN=function(x) {
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
                chrSegmented <- apply(chrSegmented, MARGIN=2L, FUN=rep,
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

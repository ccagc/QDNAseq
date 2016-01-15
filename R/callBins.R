#########################################################################/**
# @RdocFunction callBins
#
# @alias callBins,QDNAseqCopyNumbers-method
#
# @title "Call aberrations from segmented copy number data"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{object}{An object of class QDNAseqCopyNumbers}
#     \item{organism}{Either \dQuote{human} or \dQuote{other}, see manual page
#         for @see "CGHcall::CGHcall" for more details. This is only used for
#         chromosome arm information when \dQuote{prior} is set to \dQuote{all}
#         or \dQuote{auto} (and samplesize > 20). Ignored when \code{method} is
#         not \dQuote{CGHcall}.}
#     \item{method}{Calling method to use. Options currently implemented are:
#         \dQuote{CGHcall} or \dQuote{cutoff}.}
#     \item{cutoffLoss}{When method=\dQuote{cutoff}, the (log2-transformed)
#         threshold to use for losses. Default is -0.8. (The expected log2-ratio
#         for a single copy loss present in all cells is -1.0.)}
#     \item{cutoffGain}{When method=\dQuote{cutoff}, the (log2-transformed)
#         threshold to use for gains. Default is 0.5. (The expected log2-ratio
#         for a single copy gain present in all cells is 0.6.)}
#     \item{...}{Additional arguments passed to @see "CGHcall::CGHcall".}
# }
#
# \details{
#     Chromosomal aberrations are called with \pkg{CGHcall}. It has been
#     developed for the analysis of series of cancer samples, and uses a model
#     that contains both gains and losses. If used on a single sample, or
#     especially only on a subset of chromosomes, or especially on a single
#     non-cancer sample, it may fail. 
# }
#
# \value{
#     Returns an object of class @see "QDNAseqCopyNumbers" with calling
#     results added.
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
# copyNumbersCalled <- callBins(copyNumbersSegmented)
# }
#
# @author "IS"
#
# \seealso{
#     Internally, @see "CGHcall::CGHcall" and @see "CGHcall::ExpandCGHcall" of
#     the \pkg{CGHcall} package are used.
# }
#
# @keyword manip
#*/#########################################################################
setMethod('callBins', signature=c(object='QDNAseqCopyNumbers'),
    definition=function(object, organism=c("human", "other"),
    method=c("CGHcall", "cutoff"),
    cutoffLoss=-0.8, cutoffGain=0.5,
    ...) {

    method <- match.arg(method)
    if (method == "CGHcall") {
        ## Mark van de Wiel confirms that CGHcall::CGHcall() assumes (=requires)
        ## CNs on the *log* scale. /IS (private email 'Log and non-positives'
        ## on 2013-12-18 between IS and HB).
        organism <- match.arg(organism)
        if (organism == "human") {
            chrs <- fData(object)$chromosome[binsToUse(object)]
            human <- chrs %in% c(1:22, "X", "Y", "MT")
            if (any(!human)) {
                warning(paste0("Non-human chromosome names detected:\n",
                    paste(unique(chrs[!human]), collapse=", "), ".\n",
                    "Passing 'organism=\"other\"' to CGHcall()."))
                organism <- "other"
                seg <- makeCgh(object, chromosomeReplacements="auto")
            } else {
                seg <- makeCgh(object)
            }
        } else {
            seg <- makeCgh(object, chromosomeReplacements="auto")
        }
        tryCatch({
            listcall <- CGHcall(seg, organism=organism, ...)
        }, error=function(e) {
            stop("Command CGHcall() returned the following error message:\n",
                e, "Please contact maintainer of package CGHcall: ",
                maintainer("CGHcall"), call.=FALSE)
        })
        tryCatch({
            cgh <- ExpandCGHcall(listcall, seg)
        }, error=function(e) {
            stop("Command ExpandCGHcall() returned the following error ",
                "message:\n", e,
                "Please contact maintainer of package CGHcall: ",
                maintainer("CGHcall"), call.=FALSE)
        })
        calls(object) <- calls(cgh)
        if ('probdloss' %in% assayDataElementNames(cgh)) {
            probdloss(object) <- probdloss(cgh)
        } else {
            if ('probdloss' %in% assayDataElementNames(object))
                probdloss(object) <- NULL
        }
        probloss(object) <- probloss(cgh)
        probnorm(object) <- probnorm(cgh)
        probgain(object) <- probgain(cgh)
        if ('probamp' %in% assayDataElementNames(cgh)) {
            probamp(object) <- probamp(cgh)
        } else {
            if ('probamp' %in% assayDataElementNames(object))
                probamp(object) <- NULL
        }
    } else if (method == "cutoff") {
        segmentedMatrix <- log2adhoc(assayDataElement(object, "segmented"))
        callsMatrix <- (segmentedMatrix > cutoffGain) * 1L
        callsMatrix[segmentedMatrix < cutoffLoss] <- -1L
        calls(object) <- callsMatrix
        if ("probdloss" %in% assayDataElementNames(object))
            probdloss(object) <- NULL
        probloss(object) <- (callsMatrix == -1) * 1
        probnorm(object) <- (callsMatrix == 0) * 1
        probgain(object) <- (callsMatrix == 1) * 1
        if ("probamp" %in% assayDataElementNames(object))
            probamp(object) <- NULL
    }
    object
})

# EOF

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
#     \item{cutoffs}{When method=\dQuote{cutoff}, the (log2-transformed)
#         thresholds to use for calling. Must be a numerical vector of length 4,
#         and sorted in numerical order. The first element corresponds to the
#         cutoff used for homozygous deletions, and can be NA to not distinguish
#         between hemizygous losses and homozygous deletions. The second and
#         third element represent the cutoffs used for losses and gains,
#         respectively, and cannot be NAs. The fourth element is the cutoff used
#         for amplifications, and can be NA to not distinguish between gains and
#         amplifications.}
#     \item{...}{Additional arguments passed to @see "CGHcall::CGHcall".}
# }
#
# \details{
#     By default, chromosomal aberrations are called with \pkg{CGHcall}. It has
#     been developed for the analysis of series of cancer samples, and uses a
#     model that contains both gains and losses. If used on a single sample, or
#     especially only on a subset of chromosomes, or especially on a single
#     non-cancer sample, it may fail, but method \dQuote{cutoff} can be used
#     instead.
#
#     When using method \dQuote{cutoff}, the default values assume a uniform
#     cell population and correspond to thresholds of (assuming a diploid
#     genome) 0.5, 1.5, 2.5, and 10 copies to distinguish between homozygous
#     deletions, (hemizygous) losses, normal copy number, gains, and
#     amplifications, respectively. When using with cancer samples, these values
#     might require adjustments to account for tumor cell percentage.
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
#     the \pkg{CGHcall} package are used when method=\dQuote{CGHcall}.
# }
#
# @keyword manip
#*/#########################################################################
setMethod('callBins', signature=c(object='QDNAseqCopyNumbers'),
    definition=function(object, organism=c("human", "other"),
    method=c("CGHcall", "cutoff"),
    cutoffs=log2(c(deletion=0.5, loss=1.5, gain=2.5, amplification=10) / 2),
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
        if (length(cutoffs) != 4 || !is.numeric(cutoffs))
            stop("Parameter cutoff must be a numeric vector of length 4.")
        if (any(is.na(cutoffs[2:3])))
            stop("Only first and fourth element of parameter cutoff are ",
                "allowed to be NA.")
        cutoffOrder <- sort.list(cutoffs, na.last=NA)
        if (!identical(cutoffOrder, seq_along(cutoffOrder)))
            stop("Values provided for parameter cutoff must be in numerical ",
                "order.")
        segmentedMatrix <- log2adhoc(assayDataElement(object, "segmented"))
        ## multiplication with 1L turns a logical matrix into an integer one
        ## multiplication with 1 turns a logical matrix into a numeric one
        callsMatrix <- (segmentedMatrix > cutoffs[3]) * 1L
        callsMatrix[segmentedMatrix < cutoffs[2]] <- -1L
        if (!is.na(cutoffs[1])) {
            callsMatrix[segmentedMatrix < cutoffs[1]] <- -2L
            probdloss(object) <- (callsMatrix == -2) * 1
        } else {
            if ("probdloss" %in% assayDataElementNames(object))
                probdloss(object) <- NULL
        }
        if (!is.na(cutoffs[4])) {
            callsMatrix[segmentedMatrix > cutoffs[4]] <- 2L
            probamp(object) <- (callsMatrix == 2) * 1
        } else {
            if ("probamp" %in% assayDataElementNames(object))
                probamp(object) <- NULL
        }
        calls(object) <- callsMatrix
        probloss(object) <- (callsMatrix == -1) * 1
        probnorm(object) <- (callsMatrix == 0) * 1
        probgain(object) <- (callsMatrix == 1) * 1
    }
    object
})

# EOF

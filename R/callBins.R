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
#     \item{organism}{Either ‘human’ or ‘other’, see manual page for
#         @see "CGHcall::CGHcall" for more details. This is only used for
#         chromosome arm information when ‘prior’ is set to ‘all’ or ‘auto’
#         (and samplesize > 20).}
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
    definition=function(object, organism=c("human", "other"), ...) {
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
        stop("Command ExpandCGHcall() returned the following error message:\n",
            e, "Please contact maintainer of package CGHcall: ",
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
    object
})

# EOF

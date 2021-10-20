#########################################################################/**
# @RdocFunction correctBins
#
# @alias correctBins,QDNAseqReadCounts-method
#
# @title "Correct binned read counts for GC content and mappability"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{object}{An @see "QDNAseqReadCounts" object with \code{counts} data.}
#     \item{fit}{An optional matrix of values to use for the correction. If
#         NULL (default), assay data \code{fit} from object is used. If it is
#         missing, it is generated with a call to @see "estimateCorrection".}
#     \item{method}{A character(1) string specifying the correction method.
#         \code{ratio} (default) divides \code{counts} with \code{fit}.
#         \code{median} calculates the median \code{fit}, and defines the
#         correction for bins with GC content \code{gc} and mappability
#         \code{map} as \code{median(fit) - fit(gc,map)}, which is added to
#         \code{counts}. Method \code{none} leaves \code{counts} untouched.}
#     \item{adjustIncompletes}{A boolean(1) specifying whether \code{counts} for
#         bins with uncharacterized nucleotides (N's) in their reference genome
#         sequence should be adjusted by dividing them with the percentage of
#         characterized (A, C, G, T) nucleotides. Defaults to @TRUE.}
#     \item{...}{Additional arguments passed to @see "estimateCorrection".}
# }
#
# \value{
#     Returns a @see "QDNAseqCopyNumbers" object with assay data element
#     \code{copynumber}.
# }
#
# \examples{
# data(LGG150)
# readCounts <- LGG150
# readCountsFiltered <- applyFilters(readCounts)
# readCountsFiltered <- estimateCorrection(readCountsFiltered)
# copyNumbers <- correctBins(readCountsFiltered)
# }
#
# @author "IS"
#
# \seealso{
#     Internally, @see "stats::loess" is used to fit the regression model.
# }
# @keyword manip
# @keyword loess
#*/#########################################################################

setMethod("correctBins", signature=c(object="QDNAseqReadCounts"),
    definition=function(object, fit=NULL,
    method=c("ratio", "median", "none"), adjustIncompletes=TRUE, ...) {
    counts <- assayDataElement(object, "counts")
    if (adjustIncompletes) {
        counts <- counts / fData(object)$bases * 100L
        counts[fData(object)$bases == 0] <- 0L
    }
    method <- match.arg(method)
    if (method == "none")
        fit <- matrix(1, nrow=nrow(counts), ncol=ncol(counts),
            dimnames=dimnames(counts))
    if (is.null(fit)) {
        if (! "fit" %in% assayDataElementNames(object))
            object <- estimateCorrection(object, ...)
        fit <- assayDataElement(object, "fit")
    }
    if (!is.matrix(fit))
        stop("Argument fit has to be either a matrix, ",
            "or NULL, in which case estimateCorrection() is executed first.")
    if (method == "median") {
        gc <- round(fData(object)$gc)
        mappability <- round(fData(object)$mappability)
        fit2 <- aggregate(fit, by=list(gc=gc,
            mappability=mappability), FUN=median)
        rownames(fit2) <- paste0(fit2$gc, "-", fit2$mappability)
        corrected <- matrix(NA_real_, nrow=nrow(counts), ncol=ncol(counts),
            dimnames=dimnames(counts))
        for (s in sampleNames(object)) {
            correction <- median(fit2[, s], na.rm=TRUE) - fit2[, s]
            names(correction) <- rownames(fit2)
            corrected[, s] <- counts[, s] + correction[
                paste0(gc, "-", mappability)]
            corrected[, s] <- corrected[, s] - min(corrected[, s], na.rm=TRUE)
        }
    } else {
        corrected <- counts / fit
        corrected[fit <= 0] <- 0
    }
    if (anyNA(corrected[binsToUse(object), ])) {
      old <- binsToUse(object)
      missings <- is.na(rowMeans(corrected))
      message("Note: Filtering out additional ", sum(old & missings),
        " bins due to missing values.")
      binsToUse(object) <- old & !missings
    }
    new("QDNAseqCopyNumbers", bins=featureData(object), copynumber=corrected,
        phenodata=phenoData(object))
})

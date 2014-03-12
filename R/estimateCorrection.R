#########################################################################/**
# @RdocFunction estimateCorrection
#
# @alias estimateCorrection,QDNAseqReadCounts-method
#
# @title "Estimate correction to read counts for GC content and mappability"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{object}{An @see "QDNAseqReadCounts" object with \code{counts} data.}
#     \item{span}{For @see "stats::loess", the parameter alpha which controls
#         the degree of smoothing.}
#     \item{family}{For @see "stats::loess", if "gaussian" fitting is by
#         least-squares, and if "symmetric" a re-descending M estimator is
#         used with Tukey's biweight function.}
#     \item{adjustIncompletes}{A boolean(1) specifying whether \code{counts} for
#         bins with uncharacterized nucleotides (N's) in their reference genome
#         sequence should be adjusted by dividing them with the percentage of
#         characterized (A, C, G, T) nucleotides. Defaults to @TRUE.}
#     \item{maxIter}{An integer(1) specifying the maximum number of iterations
#         to perform, default is 1. If larger, after the first loess fit, bins
#         with median residuals larger than \code{cutoff} are removed, and the
#         fitting repeated until the list of bins to use stabilizes or after
#         \code{maxIter} iterations.}
#     \item{cutoff}{A numeric(1) specifying the number of standard deviations
#         (as estimated with @see "matrixStats::madDiff") the cutoff for removal
#         of bins with median residuals larger than the cutoff. Not used if
#         \code{maxIter=1} (default).}
#     \item{...}{Additional aguments passed to @see "stats::loess".}
# }
#
# \value{
#     Returns a @see "QDNAseqReadCounts" object with the assay data element
#     \code{fit} added.
# }
#
# \examples{
# data(LGG150)
# readCounts <- LGG150
# readCountsFiltered <- applyFilters(readCounts)
# readCountsFiltered <- estimateCorrection(readCountsFiltered)
# }
#
# @author "IS"
#
# \seealso{
#     Internally, @see "stats::loess" is used to fit the regression model.
# }
#
# @keyword manip
# @keyword loess
#*/#########################################################################

setMethod("estimateCorrection", signature=c(object="QDNAseqReadCounts"),
    definition=function(object, span=0.65, family="symmetric",
    adjustIncompletes=TRUE, maxIter=1, cutoff=4.0, ...) {
    
    counts <- assayDataElement(object, "counts")
    if (adjustIncompletes) {
        counts <- counts / fData(object)$bases * 100L
        counts[fData(object)$bases == 0] <- 0L
    }
    vmsg("Calculating correction for GC content and mappability:")
    if (length(span) == 1L)
        span <- rep(span, times=ncol(counts))
    if (length(family) == 1L)
        family <- rep(family, times=ncol(counts))
    if (length(span) != ncol(counts))
        stop("Parameter span has to be either a single value or a vector the ",
            "same length as there are samples in object.")
    if (length(family) != ncol(counts))
        stop("Parameter family has to be either a single value or ",
            "a vector the same length as there are samples in object.")
    condition <- binsToUse(object)
    used.span <- rep(NA_real_, times=ncol(counts))
    used.family <- rep(NA_character_, times=ncol(counts))
    loessFit <- matrix(nrow=nrow(counts), ncol=ncol(counts),
        dimnames=dimnames(counts))
    gc <- round(fData(object)$gc)
    mappability <- round(fData(object)$mappability)
    condition <- condition & !is.na(gc) & !is.na(mappability)
    # to interpolate
    all.combinations <- expand.grid(gc=unique(gc[!is.na(gc)]),
        mappability=unique(mappability[!is.na(mappability)]))
    rownames(all.combinations) <- paste0(all.combinations$gc, "-",
        all.combinations$mappability)
    for (i in seq_len(ncol(counts))) {
        if (is.na(span[i]) && is.na(family[i])) {
            vmsg("    Skipping sample ", sampleNames(object)[i], "...")
            loessFit[, i] <- 1
            next
        }
        vmsg("    Calculating fit for sample ", sampleNames(object)[i], "...")
        noProb <- FALSE
        try({
            corvals <- counts[, i]
            median.counts <- aggregate(counts[condition, i],
                by=list(gc=gc[condition], mappability=mappability[condition]),
                FUN=median)
            rownames(median.counts) <- paste0(median.counts$gc, "-",
                median.counts$mappability)
            l <- loess(x ~ gc * mappability,
                data=median.counts, span=span[i], family=family[i], ...)
            # fit <- l$fitted
            # names(fit) <- rownames(median.counts)
            fit <- as.vector(predict(l, all.combinations))
            names(fit) <- rownames(all.combinations)
            residual <- corvals / fit[paste0(gc, "-", mappability)] - 1
            cutoffValue <- cutoff * madDiff(residual, na.rm=TRUE)
            prevGoodBins <- condition
            goodBins <- binsToUse(object) & !is.na(residual) &
                abs(residual) <= cutoffValue
            iter <- 1
            while(!identical(goodBins, prevGoodBins) && iter < maxIter) {
                median.counts2 <- aggregate(counts[goodBins, i],
                    by=list(gc=gc[goodBins], mappability=mappability[goodBins]),
                    FUN=median)
                rownames(median.counts2) <- paste0(median.counts2$gc, "-",
                    median.counts2$mappability)
                l2 <- loess(x ~ gc * mappability,
                    data=median.counts2, span=span[i], family=family[i], ...)
                # fit2 <- l$fitted
                # names(fit2) <- rownames(median.counts)
                fit2 <- as.vector(predict(l2, all.combinations))
                names(fit2) <- rownames(all.combinations)
                fit[!is.na(fit2)] <- fit2[!is.na(fit2)]
                residual <- corvals / fit[paste0(gc, "-", mappability)] - 1
                prevGoodBins <- goodBins
                goodBins <- binsToUse(object) & !is.na(residual) &
                    abs(residual) <= cutoffValue
                iter <- iter + 1
            }
            loessFit[, i] <- fit[paste0(gc, "-", mappability)]
            used.span[i] <- span[i]
            used.family[i] <- family[i]
            noProb <- TRUE
        }, silent=TRUE)
        if (!noProb)
            loessFit[, i] <- 1
    }
    object$loess.span <- used.span
    object$loess.family <- used.family
    assayDataElement(object, "fit") <- loessFit
    vmsg("Done.")
    object
})
    
# EOF

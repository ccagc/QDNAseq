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
#     \item{variables}{A character vector specifying which variables to include
#         in the correction. Can be \code{c("gc", "mappability")} (the default),
#         \code{"gc"}, or \code{"mappability"}.}
#     \item{...}{Additional arguments passed to @see "stats::loess".}
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
    adjustIncompletes=TRUE, maxIter=1, cutoff=4.0,
    variables=c("gc", "mappability"), ...) {

    counts <- assayDataElement(object, "counts")
    if (adjustIncompletes) {
        counts <- counts / fData(object)$bases * 100L
        counts[fData(object)$bases == 0] <- 0L
    }
    variables <- match.arg(variables, several.ok=TRUE)
    descriptions <- c(gc="GC content", mappability="mappability")
    vmsg("Calculating correction for ",
        paste(descriptions[variables], collapse=" and "))
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
    gc <- round(fData(object)$gc)
    mappability <- round(fData(object)$mappability)
    condition <- condition & !is.na(gc) & !is.na(mappability)
    # to interpolate
    all.combinations <- expand.grid(gc=unique(gc[!is.na(gc)]),
        mappability=unique(mappability[!is.na(mappability)]))
    rownames(all.combinations) <- paste0(all.combinations$gc, "-",
        all.combinations$mappability)
    calculateFits <-  function(i, ...) {
        if (is.na(span[i]) && is.na(family[i])) {
            vmsg("    Skipping sample ", sampleNames(object)[i], "...")
            loessFit <- rep(1, times=nrow(counts))
            attr(loessFit, "used.span") <- NA_real_
            attr(loessFit, "used.family") <- NA_character_
            return(loessFit)
        }
        vmsg("    Calculating fit for sample ", sampleNames(object)[i],
          " (", i, " of ", ncol(counts), ") ...")
        try({
            corvals <- counts[, i]
            median.counts <- aggregate(counts[condition, i],
                by=list(gc=gc[condition], mappability=mappability[condition]),
                FUN=median)
            rownames(median.counts) <- paste0(median.counts$gc, "-",
                median.counts$mappability)
            l <- loess(formula(paste("x ~", paste(variables, collapse=" * "))),
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
            loessFit <- fit[paste0(gc, "-", mappability)]
            attr(loessFit, "used.span") <- span[i]
            attr(loessFit, "used.family") <- family[i]
            return(loessFit)
        }, silent=TRUE)
        loessFit <- rep(1, times=nrow(counts))
        attr(loessFit, "used.span") <- NA_real_
        attr(loessFit, "used.family") <- NA_character_
        return(loessFit)
    }
    fits <- flapply(seq_len(ncol(counts)), FUN=calculateFits, ...)
    loessFit <- do.call(cbind, fits)
    dimnames(loessFit) <- dimnames(counts)
    object$loess.span <- unlist(lapply(fits, FUN=attr, which="used.span"))
    object$loess.family <- unlist(lapply(fits, FUN=attr, which="used.family"))
    assayDataElement(object, "fit") <- loessFit
    vmsg("Done.")
    object
})

# EOF

#########################################################################/**
# @RdocFunction highlightFilters
#
# @alias highlightFilters,QDNAseqSignals-method
#
# @title "Highlights data points in a plotted profile to evaluate filtering"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{object}{A @see "QDNAseqCopyNumbers" object.}
#     \item{col}{The color used for highlighting.}
#     \item{residual}{Either a @logical specifying whether to filter based on
#         loess residuals of the calibration set, or if a @numeric, the cutoff
#         as number of standard deviations estimated with
#         @see "matrixStats::madDiff" to use for. Default is @TRUE, which
#         corresponds to 4.0 standard deviations.}
#     \item{blacklist}{Either a @logical specifying whether to filter based on
#         overlap with blacklisted regions, or if numeric, the maximum
#         percentage of overlap allowed. Default is @TRUE, which corresponds to
#         no overlap allowd (i.e. value of 0).}
#     \item{mappability}{A @numeric in \eqn{[0,100]} to specify filtering out
#         bins with mappabilities lower than the number specified. NA (default)
#         or @FALSE will not filter based on mappability.}
#     \item{bases}{A @numeric specifying the minimum percentage of characterized
#         bases (not Ns) in the reference genome sequence. NA (default) or
#         @FALSE will not filted based on uncharacterized bases.}
#     \item{type}{When specifying multiple filters (\code{residual},
#         \code{blacklist}, \code{mappability}, \code{bases}), whether to
#         highlight their \code{union} (default) or \code{intersection}.}
#     \item{...}{Further arguments to @see "graphics::points".}
# }
#
# \examples{
# data(LGG150)
# readCounts <- LGG150
# plot(readCounts)
# highlightFilters(readCounts, residual=TRUE, blacklist=TRUE)
# }
#
# @author "IS"
#
# @keyword aplot
#*/#########################################################################

setMethod("highlightFilters", signature=c(object="QDNAseqSignals"),
    definition=function(object, col="red", residual=NA, blacklist=NA,
    mappability=NA, bases=NA, type=c("union", "intersection"), ...) {

    condition <- rep(TRUE, times=nrow(object))
    type <- match.arg(type)
    if (type == "intersection") {
        if (!is.na(residual)) {
            if (is.numeric(residual)) {
                condition <- condition & (is.na(fData(object)$residual) |
                    abs(fData(object)$residual) > residual)
            } else if (residual) {
                condition <- condition & is.na(fData(object)$residual)
            }
        }
        if (!is.na(blacklist)) {
            if (is.numeric(blacklist)) {
                condition <- condition & fData(object)$blacklist > blacklist
            } else if (blacklist) {
                condition <- condition & fData(object)$blacklist != 0
            }
        }
        if (!is.na(mappability) && mappability != FALSE)
            condition <- condition & fData(object)$mappability < mappability
        if (!is.na(bases) && bases != FALSE)
            condition <- condition & fData(object)$bases < bases
    } else {
        if (!is.na(residual)) {
            if (is.numeric(residual)) {
                condition <- condition & !is.na(fData(object)$residual) &
                    abs(fData(object)$residual) <= residual
            } else if (residual) {
                condition <- condition & !is.na(fData(object)$residual)
            }
        }
        if (!is.na(blacklist)) {
            if (is.numeric(blacklist)) {
                condition <- condition & fData(object)$blacklist <= blacklist
            } else if (blacklist) {
                condition <- condition & fData(object)$blacklist == 0
            }
        }
        if (!is.na(mappability) && mappability != FALSE)
            condition <- condition & fData(object)$mappability >= mappability
        if (!is.na(bases) && bases != FALSE)
            condition <- condition & fData(object)$bases >= bases
        condition <- !condition
    }

    i <- 1
    chrs <- unique(fData(object)$chromosome[binsToUse(object)])
    condition <- condition & fData(object)$bases > 0 &
        fData(object)$chromosome %in% chrs

    all.chrom <- chromosomes(object)
    chrom <- all.chrom[condition]
    uni.chrom <- unique(chrom)
    if (!getOption("QDNAseq::plotScale")) {
        index <- 1:nrow(object)
        indexPos <- rep(NA_integer_, times=nrow(object))
        indexPos[binsToUse(object)] <- 1:sum(binsToUse(object))
        f <- approxfun(index, indexPos)
        pos <- f(index[condition])
    } else {
        all.chrom.lengths <- aggregate(bpend(object),
            by=list(chromosome=all.chrom), FUN=max)
        chrom.lengths <- all.chrom.lengths$x
        names(chrom.lengths) <- all.chrom.lengths$chromosome

        pos <- as.numeric(bpstart(object)[condition])
        chrom.lengths <- chrom.lengths[uni.chrom]
        chrom.num <- as.integer(factor(chrom, levels=uni.chrom, ordered=TRUE))
        uni.chrom.num <- unique(chrom.num)
        for (j in seq_along(uni.chrom)) {
            pos[chrom.num > uni.chrom.num[j]] <-
                pos[chrom.num > uni.chrom.num[j]] +
                chrom.lengths[uni.chrom[j]]
        }
    }
    if (class(object) == "QDNAseqReadCounts") {
        copynumber <- assayDataElement(object, "counts")[condition, ,
            drop=FALSE]
    } else {
        copynumber <- assayDataElement(object, "copynumber")[condition, ,
            drop=FALSE]
    }
    if (getOption("QDNAseq::plotLogTransform"))
        copynumber <- log2adhoc(copynumber)
    if (ncol(object) > 1L)
        vmsg("Multiple samples present in input, only using first sample: ",
            sampleNames(object)[1L])
    pointcol <- col
    ylim <- par("usr")[3:4]
    points(pos, copynumber[, i], cex=0.1, col=pointcol)
    amps <- copynumber[, i]
    amps[amps <= ylim[2]] <- NA_real_
    amps[!is.na(amps)] <- ylim[2] + 0.01 * (ylim[2]-ylim[1])
    dels <- copynumber[, i]
    dels[dels >= ylim[1]] <- NA_real_
    dels[!is.na(dels)] <- ylim[1] - 0.01 * (ylim[2]-ylim[1])
    par(xpd=TRUE)
    points(pos, amps, pch=24, col=pointcol, bg=pointcol, cex=0.5)
    points(pos, dels, pch=25, col=pointcol, bg=pointcol, cex=0.5)
    par(xpd=FALSE)

    num <- sum(condition)
    vmsg("Highlighted ", format(num, big.mark=","), " bins.")
    return(invisible(num))
})

# EOF

#########################################################################/**
# @RdocFunction applyFilters
#
# @alias applyFilters,QDNAseqReadCounts-method
#
# @title "Adjusts the filtering on which bins are used"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{object}{A @see "QDNAseqReadCounts" object.}
#     \item{residual}{Either a @logical specifying whether to filter based on
#         loess residuals of the calibration set, or if a @numeric, the cutoff
#         as number of standard deviations estimated with
#         @see "matrixStats::madDiff" to use for. Default is @TRUE, which
#         corresponds to 4.0 standard deviations.}
#     \item{blacklist}{Either a @logical specifying whether to filter based on
#         overlap with blacklisted regions, or if numeric, the maximum
#         percentage of overlap allowed. Default is @TRUE, which corresponds to
#         no overlap allowed (i.e. value of 0).}
#     \item{mappability}{A @numeric in \eqn{[0,100]} to specify filtering out
#         bins with mappabilities lower than the number specified. NA (default)
#         or @FALSE will not filter based on mappability.}
#     \item{bases}{A @numeric specifying the minimum percentage of characterized
#         bases (not Ns) in the reference genome sequence. NA (default) or
#         @FALSE will not filter based on uncharacterized bases.}
#     \item{chromosomes}{A @character vector specifying which chromosomes
#         to filter out. Defaults to the sex chromosomes and mitochondria,
#         i.e. \code{c("X", "Y", "MT")}.}
#     \item{verbose}{If @TRUE, verbose messages are produced.}
# }
#
# \value{
#     Returns a @see "QDNAseqReadCounts" object with updated filtering.
# }
#
# \examples{
# data(LGG150)
# readCounts <- LGG150
# readCountsFiltered <- applyFilters(readCounts)
# }
#
# @author "IS"
#
# @keyword manip
#*/#########################################################################
setMethod('applyFilters', signature=c(object='QDNAseqReadCounts'),
    definition=function(object, residual=TRUE, blacklist=TRUE, mappability=NA,
    bases=NA, chromosomes=c("X", "Y", "MT"),
    verbose=getOption("QDNAseq::verbose", TRUE)) {

    oopts <- options("QDNAseq::verbose"=verbose)
    on.exit(options(oopts))

    condition <- rep(TRUE, times=nrow(object))
    msg <- c('total bins'=sum(condition))
    condition <- !fData(object)$chromosome %in% chromosomes
    msg <- c(msg, 'of which in selected chromosomes'=sum(condition))
    condition <- condition & !is.na(fData(object)$gc)
    msg <- c(msg, 'of which with reference sequence'=sum(condition))

    if (!is.na(residual) && !"residual" %in% colnames(fData(object))) {
        if (is.numeric(residual) || residual) {
            warning("Residuals missing from bin annotations, filter not used.")
        }
        residual <- NA
    }
    if (!is.na(blacklist) && !"blacklist" %in% colnames(fData(object))) {
        if (is.numeric(blacklist) || blacklist) {
            warning("Blacklisted regions missing from bin annotations, ",
                "filter not used.")
        }
        blacklist <- NA
    }
    if (!is.na(mappability) && !"mappability" %in% colnames(fData(object))) {
        warning("Mappabilities missing from bin annotations, filter not used.")
        mappability <- NA
    }
    if (!is.na(bases) && !"bases" %in% colnames(fData(object))) {
        warning("Percentages of characterized bases missing from bin ",
            "annotations, filter not used.")
        bases <- NA
    }

    if (!is.na(residual)) {
        residuals <- fData(object)$residual
        cutoff <- residual * madDiff(residuals, na.rm=TRUE)
        residualsMissing <- aggregate(residuals,
            by=list(chromosome=fData(object)$chromosome),
            FUN=function(x) all(is.na(x)))
        chromosomesWithResidualsMissing <-
            residualsMissing$chromosome[residualsMissing$x]
        chromosomesToInclude <-
            setdiff(chromosomesWithResidualsMissing, chromosomes)
        if (length(chromosomesToInclude) > 0) {
            message("Note: Residual filter missing for chromosomes: ",
                paste(chromosomesToInclude, collapse=", "))
            residuals[fData(object)$chromosome %in% chromosomesToInclude] <- 0
        }
        if (is.numeric(residual)) {
            condition <- condition & !is.na(residuals) &
                abs(residuals) <= cutoff
        } else if (residual) {
            condition <- condition & !is.na(residuals)
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
    msg <- c(msg, 'final bins'=sum(condition))

    binsToUse(object) <- condition
    ## Assigning to an 'eSet' object will drop dimnames, so we can equally
    ## well use useNames=FALSE here to minimize memory and speed overheads
    object$used.reads <- colSums2(assayDataElement(object, "counts"), rows=condition, useNames=FALSE)
    object$expected.variance <- expectedVariance(object)
    vmsg(paste(format(msg, big.mark=','), names(msg),
        sep='\t', collapse='\n'))
    object
})

# EOF

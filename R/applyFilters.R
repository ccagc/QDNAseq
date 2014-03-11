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
#  @get "title".
# }
#
# \arguments{
#   \item{object}{A @see "QDNAseqReadCounts" object.}
#   \item{residual}{Either a @logical specifying whether to filter based on
#     loess residuals of the calibration set, or if a @numeric, the cutoff as
#     number of standard deviations estimated with @see "matrixStats::madDiff"
#     to use for. Default is @TRUE, which corresponds to 4.0 standard
#     deviations.}
#   \item{blacklist}{Either a @logical specifying whether to filter based on
#     overlap with blacklisted regions, or if numeric, the maximum percentage of
#     overlap allowed. Default is @TRUE, which corresponds to no overlap allowd
#     (i.e. value of 0).}
#   \item{mappability}{A @numeric in \eqn{[0,100]} to specify filtering out bins
#     with mappabilities lower than the number specified. NA (default) or
#     @FALSE will not filter based on mappability.}
#   \item{bases}{A @numeric specifying the minimum percentage of characterized
#     bases (not Ns) in the reference genome sequence. NA (default) or
#     @FALSE will not filted based on uncharacterized bases.}
#   \item{filterAllosomes}{A @logical specifying whether to filter out the sex
#     chromosomes. Default is @TRUE.}
# }
#
# \value{
#   Returns a @see "QDNAseqReadCounts" object with updated filtering.
# }
#
# \examples{
#   data(LGG150)
#   readCounts <- LGG150
#   readCountsFiltered <- applyFilters(readCounts)
# }
#
# @author "IS"
#
# @keyword IO
#*/#########################################################################
setMethod('applyFilters', signature=c(object='QDNAseqReadCounts'),
  definition=function(object, residual=TRUE, blacklist=TRUE, mappability=NA,
  bases=NA, filterAllosomes=TRUE) {

  condition <- rep(TRUE, times=nrow(object))
  msg <- c('total bins'=sum(condition))
  if (filterAllosomes) {
    condition <- fData(object)$chromosome %in% as.character(1:22)
    msg <- c(msg, 'autosomal bins'=sum(condition))
  }
  condition <- condition & !is.na(fData(object)$gc)
  msg <- c(msg, 'bins with reference sequence'=sum(condition))

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
  msg <- c(msg, 'final bins'=sum(condition))

  binsToUse(object) <- condition
  object$used.reads <- colSums(assayDataElement(object, "counts")[condition, ,
    drop=FALSE])
  object$expected.variance <- sum(condition) / object$used.reads
  vmsg(paste(format(msg, big.mark=','), names(msg),
    sep='\t', collapse='\n'))
  object
})

# EOF

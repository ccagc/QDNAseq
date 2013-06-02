#########################################################################/**
# @RdocFunction filterBins
#
# @title "Adjusts the filter on which bins are used"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{obj}{A qdnaseq object ...}
#   \item{blacklist}{...}
#   \item{mappability}{...}
#   \item{tgr}{...}
#   \item{bases}{...}
#   \item{allosomes}{...}
# }
#
# \value{
#   Returns a named @list with element 'phenodata' ...
# }
#
# @author "IS"
#
# \seealso{
#   ...
# }
#
# @keyword IO
#*/#########################################################################
filterBins <- function(obj, blacklist=0, mappability=50, tgr=3, bases=100, allosomes=TRUE) {
  condition <- condition <- rep(TRUE, times=nrow(obj))
  condition <- condition & fData(obj)$bases >= bases
  condition <- condition & fData(obj)$blacklist <= blacklist
  condition <- condition & fData(obj)$mappability >= mappability
  condition <- condition & fData(obj)$tgr <= tgr*sd(fData(obj)$tgr, na.rm=TRUE)
  if (allosomes)
    condition <- condition & fData(obj)$chromosome %in% as.character(1:22)
  fData(obj)$filter <- condition
  cat('Total bins:    ', format(nrow(obj), trim=TRUE, big.mark=','), '\n', sep='')
  cat('Bins filtered: ', format(sum(!condition), trim=TRUE, big.mark=','), '\n', sep='')
  cat('Final bins:    ', format(sum(condition), trim=TRUE, big.mark=','), '\n', sep='')
  obj
}

# EOF

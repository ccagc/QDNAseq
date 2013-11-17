#########################################################################/**
# @RdocFunction applyFilters
#
# @alias applyFilters,QDNAseqReadCounts-method
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
#   \item{object}{A QDNAseqReadCounts object ...}
#   \item{residual}{...}
#   \item{blacklist}{...}
#   \item{mappability}{...}
#   \item{bases}{...}
#   \item{filterAllosomes}{...}
#   \item{force}{...}
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
setMethod('applyFilters', signature=c(object='QDNAseqReadCounts'),
  definition=function(object, residual=TRUE, blacklist=TRUE, mappability=NA,
  bases=NA, filterAllosomes=TRUE, force=FALSE) {
  if (!force && 'segmented' %in% assayDataElementNames(object))
    stop('Data has already been segmented. Changing the filters will ',
      'remove segmentation (and possible calling) results. Please specify ',
      'force=TRUE, if you want this.')
  if ('segmented' %in% assayDataElementNames(object))
    assayDataElement(object, 'segmented') <- NULL
  if ('calls' %in% assayDataElementNames(object)) {
    assayDataElement(object, 'calls') <- NULL
    assayDataElement(object, 'probloss') <- NULL
    assayDataElement(object, 'probnorm') <- NULL
    assayDataElement(object, 'probgain') <- NULL
    if ('probdloss' %in% assayDataElementNames(object))
      assayDataElement(object, 'probdloss') <- NULL
    if ('probamp' %in% assayDataElementNames(object))
      assayDataElement(object, 'probamp') <- NULL
  }

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
  vmsg(paste(format(msg, big.mark=','), names(msg),
    sep='\t', collapse='\n'))
  object
})

# EOF

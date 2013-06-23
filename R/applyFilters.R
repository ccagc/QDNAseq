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
#   \item{mappability}{...}
#   \item{blacklist}{...}
#   \item{residual}{...}
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
  definition=function(object, mappability=50, blacklist=0, residual=2,
  bases=100, filterAllosomes=TRUE, force=FALSE) {
  if ('segmented' %in% assayDataElementNames(object) & !force)
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
  condition <- condition <- rep(TRUE, times=nrow(object))
  if (!is.na(bases))
    condition <- condition & fData(object)$bases >= bases
  if (!is.na(blacklist))
    condition <- condition & fData(object)$blacklist <= blacklist
  if (!is.na(mappability))
    condition <- condition & fData(object)$mappability >= mappability
  if (!is.na(residual))
    condition <- condition & abs(fData(object)$residual) <=
      residual*sd(fData(object)$residual, na.rm=TRUE)
  if (filterAllosomes) {
    condition2 <- fData(object)$chromosome %in% as.character(1:22)
    condition <- condition & condition2
    message('Flagging allosomes for filtering. ',
      'Statistics for autosomes 1-22:')
  } else {
    condition2 <- rep(TRUE, times=nrow(object))
    message('Statistics for all chromosomes:')
  }
  fData(object)$filter <- condition
  message(paste(c('Total bins:', 'Bins filtered:', 'Final bins:'),
    format(c(nrow(object[condition2,]),
    sum(!condition[condition2], na.rm=TRUE),
    sum(condition[condition2], na.rm=TRUE)),
    big.mark=','), sep='\t', collapse='\n'))
  object
})

# EOF

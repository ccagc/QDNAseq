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
#   \item{blacklist}{...}
#   \item{mappability}{...}
#   \item{tgr}{...}
#   \item{bases}{...}
#   \item{allosomes}{...}
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
setMethod('applyFilters', signature=c(object='QDNAseqReadCounts'), definition=function(object,
  blacklist=0, mappability=50, tgr=2, bases=100, allosomes=TRUE,
  force=FALSE) {
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
  condition <- condition & fData(object)$bases >= bases
  condition <- condition & fData(object)$blacklist <= blacklist
  condition <- condition & fData(object)$mappability >= mappability
  condition <- condition & abs(fData(object)$tgr) <= tgr*sd(fData(object)$tgr,
    na.rm=TRUE)
  if (allosomes) {
    condition2 <- fData(object)$chromosome %in% as.character(1:22)
    condition <- condition & condition2
    message('Excluding autosomes.')
  } else {
    condition2 <- rep(TRUE, times=nrow(object))
  }
  fData(object)$filter <- condition
  message('Total bins:\t', format(nrow(object[condition2,]),
    trim=TRUE, big.mark=','))
  message('Bins filtered:\t', format(sum(!condition[condition2], na.rm=TRUE),
    trim=TRUE, big.mark=','))
  message('Final bins:\t', format(sum(condition[condition2], na.rm=TRUE),
    trim=TRUE, big.mark=','))
  object
})

# EOF

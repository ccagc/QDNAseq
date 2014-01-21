#########################################################################/**
# @RdocFunction makeCgh
#
# @alias makeCgh,QDNAseqReadCounts-method
#
# @title "Constructs a 'cghRaw', 'cghSeg', or 'cghCal' object"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{object}{...}
# }
#
# \value{
#   Returns a @see "CGHbase::cghRaw" object.
# }
#
# @author "IS"
#*/#########################################################################
setMethod('makeCgh', signature=c(object='QDNAseqReadCounts'),
  definition=function(object) {
  object <- object[binsToUse(object),]
  fData(object)$chromosome <- chromosomes(object)
  colnames(fData(object))[colnames(fData(object)) == 'chromosome'] <-
    'Chromosome'
  colnames(fData(object))[colnames(fData(object)) == 'start'] <-
    'Start'
  colnames(fData(object))[colnames(fData(object)) == 'end'] <-
    'End'
  if ('counts' %in% assayDataElementNames(object))
    assayDataElement(object, 'counts') <- NULL
  if ('corrected' %in% assayDataElementNames(object))
    assayDataElement(object, 'corrected') <- NULL
  if ('residuals' %in% assayDataElementNames(object))
    assayDataElement(object, 'residuals') <- NULL
  if ('calls' %in% assayDataElementNames(object)) {
    className <- 'cghCall'
  } else if ('segmented' %in% assayDataElementNames(object)) {
    className <- 'cghSeg'
  } else {
    className <- 'cghRaw'
  }
  cgh <- new(className, assayData=assayData(object),
    featureData=featureData(object), phenoData=phenoData(object))
  cgh
})

# EOF

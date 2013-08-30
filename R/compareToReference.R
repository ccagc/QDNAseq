#########################################################################/**
# @RdocFunction compareToReference
#
# @alias compareToReference,QDNAseqReadCounts,numeric-method
#
# @title "Divide binned read counts with those of reference samples"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{object}{...}
#   \item{references}{...}
#   \item{force}{...}
# }
#
# \value{
#   Returns a named @list containing elements ...
# }
#
# @author "IS"
#
# \seealso{
#   Internally, ...
# }
#
#*/#########################################################################
setMethod('compareToReference', signature=c(object='QDNAseqReadCounts',
  references='numeric'), definition=function(object, references,
  force=FALSE) {

  if (!force && 'copynumber' %in% assayDataElementNames(object))
    stop('Data has already been normalized. Comparing to reference will ',
      'remove normalization (and possible segmentation and calling) ',
      'results. Please specify force=TRUE, if you want this.')
  if ('copynumber' %in% assayDataElementNames(object))
    assayDataElement(object, 'copynumber') <- NULL
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

  if (length(references) != ncol(object))
    stop('Parameter references must be a vector of equal length as there ',
      'are samples in object.')
  for (i in seq_along(references)) {
    if (!is.na(references[i]) && references[i] != FALSE) {
      assayDataElement(object, 'counts')[, i] <-
        (assayDataElement(object, 'counts')[, i]+1) /
        (assayDataElement(object, 'counts')[, references[i]]+1) - 1
      if ('corrected' %in% assayDataElementNames(object))
        assayDataElement(object, 'corrected')[, i] <-
          (assayDataElement(object, 'corrected')[, i]+1) /
          (assayDataElement(object, 'corrected')[, references[i]]+1) - 1
      sampleNames(object)[i] <- paste(sampleNames(object)[i], ' vs. ',
        sampleNames(object)[references[i]], sep='')
    }
  }
  toremove <- which(!references)
  if (length(toremove) > 0L)
    object <- object[, -toremove]
  object
})

# EOF

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
  references='numeric'), definition=function(object, references) {
  if (length(references) != ncol(object))
    stop('Parameter references must be a vector of equal length as there ',
      'are samples in object.')
  for (i in seq_along(references)) {
    if (!is.na(references[i]) && references[i]!=FALSE) {
      assayDataElement(object, 'counts')[, i] <- (assayDataElement(object,
        'counts')[, i]+1) / (assayDataElement(object,
        'counts')[, references[i]]+1) - 1
      if ('corrected' %in% assayDataElementNames(object))
        assayDataElement(object, 'corrected')[, i] <- (assayDataElement(object,
          'corrected')[, i]+1) / (assayDataElement(object,
          'corrected')[, references[i]]+1)
      sampleNames(object)[i] <- paste(sampleNames(object)[i], ' vs. ',
        sampleNames(object)[references[i]], sep='')
    }
  }
  toremove <- which(references == FALSE)
  if (length(toremove)>0)
    object <- object[, -toremove]
  object
})

# EOF

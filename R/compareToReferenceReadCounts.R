#########################################################################/**
# @RdocFunction compareToReferenceReadCounts
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
#   \item{obj}{...}
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
compareToReferenceReadCounts <- function(obj, references) {
  if (length(references) != ncol(obj))
    stop('Parameter references must be a vector of equal length as there ',
      'are samples in obj.')
  for (i in seq_along(references)) {
    if (!is.na(references[i]) && references[i]!=FALSE) {
      assayDataElement(obj, 'counts')[,i] <- (assayDataElement(obj,
        'counts')[,i]+1) / (assayDataElement(obj, 'counts')[,references[i]]+1)
        - 1
      if ('corrected' %in% assayDataElementNames(obj))
        assayDataElement(obj, 'corrected')[,i] <- (assayDataElement(obj,
          'corrected')[,i]+1) / (assayDataElement(obj,
          'corrected')[,references[i]]+1)
      sampleNames(obj)[i] <- paste(sampleNames(obj)[i], ' vs. ',
        sampleNames(obj)[references[i]], sep='')
    }
  }
  toremove <- which(references==FALSE)
  if (length(toremove)>0)
    obj <- obj[,-toremove]
  obj
}

# EOF

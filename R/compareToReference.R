#########################################################################/**
# @RdocFunction compareToReference
#
# @alias compareToReference,QDNAseqCopyNumbers,numeric-method
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
setMethod("compareToReference", signature=c(object="QDNAseqCopyNumbers",
  references="numeric"), definition=function(object, references, force=FALSE) {

  if (!force && "segmented" %in% assayDataElementNames(object))
    stop("Data has already been segmented. Comparing to reference will ",
      "remove segmentation (and possible calling) ",
      "results. Please specify force=TRUE, if you want this.")
  if ("segmented" %in% assayDataElementNames(object))
    assayDataElement(object, "segmented") <- NULL
  if ("calls" %in% assayDataElementNames(object)) {
    assayDataElement(object, "calls") <- NULL
    assayDataElement(object, "probloss") <- NULL
    assayDataElement(object, "probnorm") <- NULL
    assayDataElement(object, "probgain") <- NULL
    if ("probdloss" %in% assayDataElementNames(object))
      assayDataElement(object, "probdloss") <- NULL
    if ("probamp" %in% assayDataElementNames(object))
      assayDataElement(object, "probamp") <- NULL
  }

  if (length(references) != ncol(object))
    stop("Parameter references must be a vector of equal length as there ",
      "are samples in object.")
  for (i in seq_along(references)) {
    if (!is.na(references[i]) && references[i] != FALSE) {
      test <- assayDataElement(object, "copynumber")[, i]
      reference <- assayDataElement(object, "copynumber")[, references[i]]
      ratio <- test / reference
      ratio[reference == 0] <- 0
      assayDataElement(object, "copynumber")[, i] <- ratio
      sampleNames(object)[i] <- paste(sampleNames(object)[i], " vs. ",
        sampleNames(object)[references[i]], sep="")
    }
  }
  toremove <- which(!references)
  if (length(toremove) > 0L)
    object <- object[, -toremove]
  object
})

# EOF

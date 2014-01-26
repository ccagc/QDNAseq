#########################################################################/**
# @RdocFunction makeCgh
#
# @alias makeCgh,QDNAseqReadCounts-method
#
# @title "Constructs a 'cghRaw', 'cghSeg', or 'cghCall' object"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{object}{A @see "QDNAseqReadCounts" object.}
#   \item{filter}{If @TRUE, bins are filtered, otherwise not.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @see "CGHbase::cghRaw" if the object has not been segmented,
#   a @see "CGHbase::cghSeg" if it has been segmented but not called,
#   or @see "CGHbase::cghCall" if it has been called as well.
# }
#
# @author "IS"
#*/#########################################################################
setMethod('makeCgh', signature=c(object='QDNAseqReadCounts'),
  definition=function(object, filter=TRUE, ...) {

  # Decide which cgh* class to create
  names <- assayDataElementNames(object)
  if ('calls' %in% names) {
    className <- 'cghCall'
  } else if ('segmented' %in% names) {
    className <- 'cghSeg'
  } else if ('copynumber' %in% names) {
    className <- 'cghRaw'
  } else {
    stop("Cannot create a CGHbase::cgh* object without at least one of the assay data elements being 'calls', 'segmented' or 'copynumber': ", paste(sQuote(names), collapse=", "))
  }

  # Filter bins?
  if (filter) {
    keep <- binsToUse(object)
    object <- object[keep,]
    keep <- NULL # Not needed anymore
  }

  # Coerce chromosomes to integer indices
  fData(object)$chromosome <- chromosomes(object)

  # Update column names
  names <- colnames(fData(object))
  names[names == 'chromosome'] <- 'Chromosome'
  names[names == 'start'] <- 'Start'
  names[names == 'end'] <- 'End'
  colnames(fData(object)) <- names

  # Remove certain assay data elements
  names <- c('counts', 'corrected', 'residuals')
  names <- intersect(names, assayDataElementNames(object))
  for (name in names) assayDataElement(object, name) <- NULL

  # Instantiate choose cgh* object
  cgh <- new(className, assayData=assayData(object),
             featureData=featureData(object), phenoData=phenoData(object))

  cgh
})

# EOF

#########################################################################/**
# @RdocFunction makeCghRaw
#
# @title "Constructs a 'cghRaw' object"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{obj}{...}
# }
#
# \value{
#   Returns a @see "CGHbase::cghRaw" object.
# }
#
# @author "IS"
#
# \seealso{
#   Internally, @see "CGHbase::make_cghRaw" is used.
# }
#*/#########################################################################
makeCghRaw <- function(obj) {
  if ('filter' %in% colnames(fData(obj))) {
    condition <- fData(obj)$filter
  } else {
    condition <- rep(TRUE, times=nrow(obj))
  }
  cgh <- make_cghRaw(data.frame(bin=featureNames(obj)[condition], fData(obj)[condition, c('chromosome', 'start', 'end'),], assayDataElement(obj, 'copynumber')[condition, , drop=FALSE], check.names=FALSE, stringsAsFactors=FALSE))
  pData(cgh) <- pData(obj)
  fData(cgh)$bases <- bins[condition, 'bases']
  fData(cgh)$gc <- bins[condition, 'gc']
  fData(cgh)$mappability <- bins[condition, 'mappability']
  fData(cgh)$blacklist <- bins[condition, 'blacklist']
  fData(cgh)$tgr <- bins[condition, 'tgr']
  cgh
}

# EOF

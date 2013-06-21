#########################################################################/**
# @RdocFunction callBins
#
# @alias callBins,QDNAseqReadCounts-method
#
# @title "Call aberrations from copy number data"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{object}{...}
#   \item{...}{...}
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
setMethod('callBins', signature=c(object='QDNAseqReadCounts'),
  definition=function(object, ...) {
  seg <- makeCgh(object)
  listcall <- CGHcall(seg, ...)
  cgh <- ExpandCGHcall(listcall, seg)
  calls(object) <- calls(cgh)
  if ('probdloss' %in% assayDataElementNames(cgh)) {
    probdloss(object) <- probdloss(cgh)
  } else {
    if ('probdloss' %in% assayDataElementNames(object))
      probdloss(object) <- NULL
  }
  probloss(object) <- probloss(cgh)
  probnorm(object) <- probnorm(cgh)
  probgain(object) <- probgain(cgh)
  if ('probamp' %in% assayDataElementNames(cgh)) {
    probamp(object) <- probamp(cgh)
  } else {
    if ('probamp' %in% assayDataElementNames(object))
      probamp(object) <- NULL
  }
  object
})

# EOF

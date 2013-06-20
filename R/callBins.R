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
  probdloss(object) <- probdloss(cgh)
  probloss(object) <- probloss(cgh)
  probnorm(object) <- probnorm(cgh)
  probgain(object) <- probgain(cgh)
  probamp(object) <- probamp(cgh)
  object
})

# EOF

#########################################################################/**
# @RdocFunction CGHcall
#
# @alias CGHcall,QDNAseqReadCounts-method
# @alias CGHcall,cghSeg-method
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
setMethod('CGHcall', signature=c(object='QDNAseqReadCounts'),
  definition=function(object, ...) {
  seg <- makeCgh(object)
  CGHcall::CGHcall(seg, ...)
})

setMethod('CGHcall', signature=c(object='cghSeg'),
  definition=function(object, ...) {
  CGHcall::CGHcall(object, ...)
})

# EOF

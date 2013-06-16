#########################################################################/**
# @RdocFunction ExpandCGHcall
#
# @alias ExpandCGHcall,list,QDNAseqReadCounts-method
# @alias ExpandCGHcall,list,cghSeg-method
#
# @title "Expands result from CGHcall to a QDNAseqReadCounts or CGHcall object"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{listcall}{...}
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
setMethod('ExpandCGHcall', signature=
  c(listcall='list', object='QDNAseqReadCounts'),
  definition=function(listcall, object, ...) {
  seg <- makeCgh(object)
  cgh <- ExpandCGHcall(listcall, seg)
  calls(object) <- calls(cgh)
  probdloss(object) <- probdloss(cgh)
  probloss(object) <- probloss(cgh)
  probnorm(object) <- probnorm(cgh)
  probgain(object) <- probgain(cgh)
  probamp(object) <- probamp(cgh)
  object
})

setMethod('ExpandCGHcall', signature=
  c(listcall='list', object='cghSeg'),
  definition=function(listcall, object, ...) {
  CGHcall::ExpandCGHcall(listcall, object, ...)
})

# EOF

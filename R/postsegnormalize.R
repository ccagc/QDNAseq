#########################################################################/**
# @RdocFunction postsegnormalize
#
# @alias postsegnormalize,QDNAseqReadCounts-method
# @alias postsegnormalize,cghSeg-method
#
# @title "Post-segmentation normalization"
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
setMethod('postsegnormalize', signature=c(object='QDNAseqReadCounts'),
  definition=function(object, ...) {
  seg <- makeCgh(object)
  psn <- postsegnormalize(seg, ...)
  copynumber(object) <- copynumber(psn)
  segmented(object) <- segmented(psn)
  object
})

setMethod('postsegnormalize', signature=c(object='cghSeg'),
  definition=function(object, ...) {
  CGHcall::postsegnormalize(object, ...)
})

# EOF

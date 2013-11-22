#########################################################################/**
# @RdocFunction callBins
#
# @alias callBins,QDNAseqReadCounts-method
#
# @title "Call aberrations from segmented copy number data"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{object}{An object of class QDNAseqReadCounts}
#   \item{...}{Additional arguments passed to @see "CGHcall::CGHcall".}
# }
#
# \value{
#   Returns an object of class QDNAseqReadCounts with segmentation results
#     added.
# }
#
# @author "IS"
#
# \seealso{
#   Internally, @see "CGHcall::CGHcall" and @see "CGHcall::ExpandCGHcall" of
#     the \pkg{DNAcopy} package are used.
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

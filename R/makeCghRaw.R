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
  if (exists('filter', obj)) {
    condition <- obj[['filter']]
  } else {
    condition <- rep(TRUE, nrow(obj[['bins']]))
  }
  cgh <- make_cghRaw(data.frame(bin=rownames(obj[['bins']][condition,]), obj[['bins']][condition, c('chromosome', 'start', 'end'),], obj[['copynumber']][condition,], check.names=FALSE, stringsAsFactors=FALSE))
  pData(cgh) <- obj[['phenodata']]
  cgh
}

# EOF

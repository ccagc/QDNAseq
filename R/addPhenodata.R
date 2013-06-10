#########################################################################/**
# @RdocFunction addPhenodata
#
# @title "Gets bin annotation data for a particular bin size"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{obj}{A named @list containing ...}
#   \item{phenofile}{A @connection of a filename specifying ...}
# }
#
# \value{
#   Returns a named @list with element 'phenodata' ...
# }
#
# @author "IS"
#
# \seealso{
#   ...
# }
#
# @keyword IO
#*/#########################################################################
addPhenodata <- function(obj, phenofile) {
  pdata <- read.table(phenofile, header=TRUE, sep='\t', as.is=TRUE,
    row.names=1L)
  pData(obj) <- cbind(pData(obj), pdata[sampleNames(obj), ])
  obj
}

# EOF

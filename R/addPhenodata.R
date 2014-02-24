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
#   \item{object}{A @see "QDNAseqReadCounts" or @see "QDNAseqCopyNumbers"
#   object.}
#   \item{phenofile}{A file name with phenotypic data for samples in
#   \code{object}.}
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
addPhenodata <- function(object, phenofile) {
  pdata <- read.table(phenofile, header=TRUE, sep='\t', as.is=TRUE,
    row.names=1L)
  pData(object) <- cbind(pData(object), pdata[sampleNames(object), ])
  object
}

# EOF

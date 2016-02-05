#########################################################################/**
# @RdocFunction addPhenodata
#
# @title "Adds phenotype data from a file to a QDNAseqReadCounts or a
#     QDNAseqCopyNumbers object"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{object}{A @see "QDNAseqReadCounts" or @see "QDNAseqCopyNumbers"
#         object.}
#     \item{phenofile}{A file name with phenotypic data for samples in
#         \code{object}.}
# }
#
# \value{
#     Returns a @see "QDNAseqReadCounts" or @see "QDNAseqCopyNumbers" object
#     with phenotype data added.
# }
#
# \examples{
# data(LGG150)
# readCounts <- LGG150
# \dontrun{
# readCounts <- addPhenodata(readCounts, "phenodata.txt")
# }
# }
# @author "IS"
#
# @keyword IO
# @keyword file
#*/#########################################################################
addPhenodata <- function(object, phenofile) {
    pdata <- read.table(phenofile, header=TRUE, sep='\t', as.is=TRUE,
        row.names=1L)
    pData(object) <- cbind(pData(object), pdata[sampleNames(object), ,
        drop=FALSE])
    object
}

# EOF

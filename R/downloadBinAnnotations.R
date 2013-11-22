#########################################################################/**
# @RdocFunction downloadBinAnnotations
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
#   \item{binSize}{A @numeric scalar specifying the width of the bins
#    in units of kbp (1000 base pairs), e.g. \code{binSize=15} corresponds
#    to 15 kbp bins.}
#   \item{genome}{A @character string specify the genome and genome version
#    to be used.}
#   \item{type}{A @character string specify the experiment type, e.g. "SR50"
#    or "PE1000".}
#   \item{force}{If @TRUE, the bin anonnation data is retrieved/calculated
#    regardless of it already exists in the cache or not.}
# }
#
# \value{
#   Returns a @see "Biobase::AnnotatedDataFrame" object.
# }
#
# @author "IS"
#
# \seealso{
#   @see "createBins".
# }
#
# @keyword IO
#*/#########################################################################
downloadBinAnnotations <- function(binSize, genome='hg19', type='SR50',
  force=FALSE) {
  vmsg('Downloading bin annotations for genome ', genome,
    ', bin size ', binSize, 'kbp, and experiment type ', type, ' ...',
    appendLF=FALSE)
  filename <- sprintf('QDNAseq.%s.%gkbp.%s.rds', genome, binSize, type)
  urlPath <- 'http://cdn.bitbucket.org/ccagc/qdnaseq/downloads'
  remotefile <- file.path(urlPath, filename, fsep='/')
  localfile <- tempfile()
  tryCatch({
    result <- downloadFile(remotefile, localfile)
  }, error=function(e) {
    vmsg(' not found.')
    stop(e)
  })
  bins <- readRDS(localfile)
  file.remove(localfile)
  vmsg()
  bins
}

# EOF
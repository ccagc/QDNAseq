#########################################################################/**
# @RdocFunction segmentDataWithWeights
#
# @title "Segments and calls total copy numbers"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{input}{...}
#   \item{weights}{Either @TRUE or a vector of weights. If @TRUE,
#     the 1000 Genomes residuals are used as weights.}
#   \item{...}{Additional arguments passed to @see "DNAcopy::segment".}
# }
#
# \value{
#   Returns ...
# }
#
# @author "IS"
#
# \seealso{
#   Internally, @see "DNAcopy::segment" of the \pkg{DNAcopy} package,
#   which implements the CBS method, is used to segment the data.
# }
#
#*/#########################################################################
segmentDataWithWeights <- function(input, weights=TRUE, ...) {
  if (length(weights)==1L & weights) {
    if (!'tgr' %in% colnames(fData(input)))
      stop('1000 Genomes residuals not found.')
    tgr <- fData(input)$tgr
    input <- input[!is.na(tgr),]
    tgr <- tgr[!is.na(tgr)]
    tgr <- abs(tgr)
    tgr[tgr==0] <- min(tgr[tgr!=0], na.rm=TRUE)
    weights <- 1/tgr
  }
  CNA.object <- CNA(genomdat=copynumber(input), chrom=chromosomes(input), maploc=bpstart(input), data.type="logratio", sampleid=sampleNames(input))
  cat("Start data segmentation .. \n")
  if (length(weights)==1L && !weights) {
    segmented <- segment(CNA.object, ...)
  } else {
    segmented <- segment(CNA.object, weights=weights, ...)
  }
  numclone <- segmented$output$num.mark
  smrat <- segmented$output$seg
  numsmrat <- cbind(smrat, numclone)
  repdata <- function(row) {
    rep(row[1L], times=row[2L])
  }
  makelist <- apply(numsmrat, MARGIN=1L, FUN=repdata)
  joined <- unlist(makelist)
  # Not needed anymore
  rm(list="makelist")
  joined <- matrix(joined, ncol=ncol(input), byrow=FALSE)
  joined <- CGHcall:::.assignNames(joined, input)
  result <- CGHcall:::.segFromRaw(input, joined)
  pData(result) <- pData(input)
  result
}

# EOF

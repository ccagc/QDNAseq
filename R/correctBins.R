#########################################################################/**
# @RdocFunction correctBins
#
# @alias correctBins,QDNAseqReadCounts-method
#
# @title "Correct binned read counts for GC content and mappability"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{object}{An @see "QDNAseqReadCounts" object with \code{counts} data.}
#   \item{fit}{...}
#   \item{type}{...}
#   \item{adjustIncompletes}{...}
#   \item{...}{Additional arguments passed to @see "estimateCorrection".}
# }
#
# \value{
#   Returns a @see "QDNAseqReadCounts" object with the assay data element
#   \code{corrected} added.
# }
#
# @author "IS"
#
# \seealso{
#   Internally, @see "stats::loess" is used to fit the regression model.
# }
#*/#########################################################################

setMethod("correctBins", signature=c(object="QDNAseqReadCounts"),
  definition=function(object, fit=NULL,
  type=c("ratio", "addition", "none"), adjustIncompletes=TRUE, ...) {
  counts <- assayDataElement(object, "counts")
  if (adjustIncompletes) {
    counts <- counts / fData(object)$bases * 100L
    counts[fData(object)$bases == 0] <- 0L
  }
  type <- match.arg(type)
  if (type == "none")
    fit <- matrix(1, nrow=nrow(counts), ncol=ncol(counts),
      dimnames=dimnames(counts))
  if (is.null(fit)) {
    if (! "fit" %in% assayDataElementNames(object))
      object <- estimateCorrection(object, ...)
    fit <- assayDataElement(object, "fit")
  }
  if (!is.matrix(fit))
    stop("Argument fit has to be either a matrix, a QDNAseqReadCounts object, ",
      "or NULL, in which case estimateCorrection() is executed first.")
  if (type == "addition") {
    gc <- round(fData(object)$gc)
    mappability <- round(fData(object)$mappability)
    fit2 <- aggregate(fit, by=list(gc=gc,
      mappability=mappability), FUN=median)
    rownames(fit2) <- paste0(fit2$gc, "-", fit2$mappability)
    corrected <- matrix(NA_real_, nrow=nrow(counts), ncol=ncol(counts),
      dimnames=dimnames(counts))
    for (s in sampleNames(object)) {
      correction <- median(fit2[, s], na.rm=TRUE) - fit2[, s]
      names(correction) <- rownames(fit2)
      corrected[, s] <- counts[, s] + correction[paste0(gc, "-", mappability)]
      corrected[, s] <- corrected[, s] - min(corrected[, s], na.rm=TRUE)
    }
  } else {
    corrected <- counts / fit
    corrected[fit < 0] <- 0
  }
  new("QDNAseqCopyNumbers", bins=featureData(object), copynumber=corrected,
    phenodata=phenoData(object))
})

# EOF

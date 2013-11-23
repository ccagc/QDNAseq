#########################################################################/**
# @RdocFunction poolRuns
#
# @alias poolRuns,QDNAseqReadCounts,character-method
#
# @title "Pools binned read counts across samples"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{object}{...}
#   \item{samples}{...}
#   \item{force}{...}
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
setMethod('poolRuns', signature=c(object='QDNAseqReadCounts',
  samples='character'), definition=function(object, samples, force=FALSE) {

  if (!force && 'copynumber' %in% assayDataElementNames(object))
    stop('Data has already been normalized. Pooling runs will ',
      'remove normalization (and possible segmentation and calling) ',
      'results. Please specify force=TRUE, if you want this.')

  if (length(samples) != ncol(object))
    stop('Parameter samples must be a vector of equal length as there are ',
      ' samples in object.')

  newsamples <- sort(unique(samples))
  if (length(newsamples) == ncol(object)) {
    sampleNames(object) <- samples
    object <- object[,order(samples)]
    return(object)
  }
  oldsamples <- sampleNames(object)

  phenodata <- pData(object)
  phenodata[, 1] <- samples
  bins <- featureData(object)
  counts <- assayDataElement(object, 'counts')
  if ('corrected' %in% assayDataElementNames(object))
    corrected <- assayDataElement(object, 'corrected')

  newphenodata <- phenodata[0, ]
  newcounts <- matrix(nrow=nrow(counts), ncol=length(newsamples),
    dimnames=list(rownames(counts), newsamples))
  if ('corrected' %in% assayDataElementNames(object))
    newcorrected <- matrix(nrow=nrow(counts), ncol=length(newsamples),
      dimnames=list(rownames(counts), newsamples))

  concatenateIfNotEqual <- function(x) {
    x <- sort(unique(x))
    paste(x, collapse=';')
  }

  for (newsample in newsamples) {
    replicates <- samples == newsample
    newcounts[, newsample] <- rowSums(counts[, replicates, drop=FALSE])
    if ('corrected' %in% assayDataElementNames(object))
      newcorrected[, newsample] <- rowMeans(corrected[, replicates,
        drop=FALSE])
    oldphenodata <- phenodata[replicates, ]
    totalReads <- sum(oldphenodata$reads)
    numericCols <- sapply(oldphenodata, is.numeric)
    oldphenodata[1, numericCols] <- apply(oldphenodata[, numericCols,
      drop=FALSE], 2, mean)
    oldphenodata[1, !numericCols] <- apply(oldphenodata[, !numericCols,
      drop=FALSE], 2, concatenateIfNotEqual)
    oldphenodata[1, 'reads'] <- totalReads
    newphenodata <- rbind(newphenodata, oldphenodata[1,])
  }
  rownames(newphenodata) <- newphenodata[, 1]

  object2 <- new('QDNAseqReadCounts', bins=bins, counts=newcounts,
    phenodata=newphenodata)
  if ('corrected' %in% assayDataElementNames(object))
    assayDataElement(object2, 'corrected') <- newcorrected
  object2
})

# EOF

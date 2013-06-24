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
  samples='character'), definition=function(object, samples) {
  phenodata <- pData(object)
  bins <- featureData(object)
  counts <- assayDataElement(object, 'counts')
  if ('corrected' %in% assayDataElementNames(object))
    corrected <- assayDataElement(object, 'corrected')
  if (length(samples) != nrow(phenodata))
    stop('Parameter samples must be a vector of equal length as there are ',
      ' samples in object.')
  newsamples <- sort(unique(samples))
  oldsamples <- phenodata$name
  if (length(newsamples) == nrow(phenodata))
    return(object)
  newphenodata <- data.frame(name=newsamples, reads=NA, loess.span=NA,
    loess.family=NA, row.names=newsamples, stringsAsFactors=FALSE)
  ## FIXME: other phenodata variables get droppped, should be kept
  newcounts <- matrix(nrow=nrow(counts), ncol=length(newsamples),
    dimnames=list(rownames(counts), newsamples))
  if ('corrected' %in% assayDataElementNames(object))
    newcorrected <- matrix(nrow=nrow(counts), ncol=length(newsamples),
      dimnames=list(rownames(counts), newsamples))
  for (newsample in newsamples) {
    replicates <- samples == newsample
    newcounts[, newsample] <- rowSums(counts[, replicates, drop=FALSE])
    if ('corrected' %in% assayDataElementNames(object))
      newcorrected[, newsample] <- rowSums(corrected[, replicates,
        drop=FALSE])
    oldphenodata <- phenodata[replicates, ]
    newphenodata[newsample, 'reads'] <- sum(oldphenodata$reads)
    if (length(unique(oldphenodata$loess.span)) == 1L)
      newphenodata[newsample, 'loess.span'] <- oldphenodata[1L, 'loess.span']
    if (length(unique(oldphenodata$loess.family)) == 1L)
      newphenodata[newsample, 'loess.family'] <- oldphenodata[1L,
        'loess.family']
  }
  object2 <- new('QDNAseqReadCounts', bins=bins, counts=newcounts,
    phenodata=newphenodata)
  if ('corrected' %in% assayDataElementNames(object))
    assayDataElement(object2, 'corrected') <- newcorrected
  object2
})

# EOF

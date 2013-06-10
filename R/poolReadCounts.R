#########################################################################/**
# @RdocFunction poolReadCounts
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
#   \item{obj}{...}
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
poolReadCounts <- function(obj, samples) {
  phenodata <- pData(obj)
  bins <- fData(obj)
  counts <- assayDataElement(obj, 'counts')
  if ('corrected' %in% assayDataElementNames(obj))
    corrected <- assayDataElement(obj, 'corrected')
  if (length(samples) != nrow(phenodata))
    stop('Parameter samples must be a vector of equal length as there are ',
      ' samples in obj.')
  newsamples <- sort(unique(samples))
  oldsamples <- phenodata$name
  if (length(newsamples) == nrow(phenodata))
    return(obj)
  newphenodata <- data.frame(name=newsamples, reads=NA, loess.span=NA,
    loess.family=NA, row.names=newsamples, stringsAsFactors=FALSE)
  ## FIXME: other phenodata variables get droppped, should be kept
  newcounts <- matrix(nrow=nrow(counts), ncol=length(newsamples),
    dimnames=list(rownames(counts), newsamples))
  if ('corrected' %in% assayDataElementNames(obj))
    newcorrected <- matrix(nrow=nrow(counts), ncol=length(newsamples),
      dimnames=list(rownames(counts), newsamples))
  for (newsample in newsamples) {
    replicates <- samples == newsample
    newcounts[, newsample] <- rowSums(counts[, replicates, drop=FALSE])
    if ('corrected' %in% assayDataElementNames(obj))
      newcorrected[, newsample] <- rowSums(corrected[, replicates,
        drop=FALSE])
    oldphenodata <- phenodata[replicates, ]
    newphenodata[newsample, 'reads'] <- sum(oldphenodata$reads)
    if (length(unique(oldphenodata$loess.span)) == 1L)
      newphenodata[newsample, 'loess.span'] <- oldphenodata[1, 'loess.span']
    if (length(unique(oldphenodata$loess.family)) == 1L)
      newphenodata[newsample, 'loess.family'] <- oldphenodata[1,
        'loess.family']
  }
  obj2 <- new('qdnaseq', bins=bins, counts=newcounts, phenodata=newphenodata)
  if ('corrected' %in% assayDataElementNames(obj))
    assayDataElement(obj2, 'corrected') <- newcorrected
  obj2
}

# EOF

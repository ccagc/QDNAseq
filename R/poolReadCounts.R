poolReadCounts <- function(obj, samples) {
  phenodata <- obj[['phenodata']]
  bins <- obj[['bins']]
  counts <- obj[['counts']]
  if ('corrected' %in% names(obj))
    corrected <- obj[['corrected']]
  if (length(samples) != nrow(phenodata))
    stop('Vector samples must be of equal length as there are samples.')
  newsamples <- sort(unique(samples))
  oldsamples <- phenodata$name
  if (length(newsamples) == nrow(phenodata))
    return(obj)
  newphenodata <- data.frame(name=newsamples, reads=NA, loess.span=NA, loess.family=NA, row.names=newsamples, stringsAsFactors=FALSE)
  # FIXME: other phenodata variables get droppped, should be kept
  newcounts <- matrix(nrow=nrow(counts), ncol=length(newsamples), dimnames=list(rownames(counts), newsamples))
  if ('corrected' %in% names(obj))
    newcorrected <- matrix(nrow=nrow(counts), ncol=length(newsamples), dimnames=list(rownames(counts), newsamples))
  for (newsample in newsamples) {
    replicates <- samples == newsample
    newcounts[,newsample] <- rowSums(counts[,replicates, drop=FALSE])
    if ('corrected' %in% names(obj))
      newcorrected[,newsample] <- rowSums(corrected[,replicates, drop=FALSE])
    oldphenodata <- obj[['phenodata']][replicates,]
    newphenodata[newsample, 'reads'] <- sum(oldphenodata$reads)
    if (length(unique(oldphenodata$loess.span)) == 1)
      newphenodata[newsample, 'loess.span'] <- oldphenodata[1, 'loess.span']
    if (length(unique(oldphenodata$loess.family)) == 1)
      newphenodata[newsample, 'loess.family'] <- oldphenodata[1, 'loess.family']
  }
  if ('corrected' %in% names(obj))
    return(list(phenodata=newphenodata, bins=bins, counts=newcounts, filter=obj[['filter']], corrected=newcorrected))
  list(phenodata=newphenodata, bins=bins, counts=newcounts, filter=obj[['filter']])
}

# EOF

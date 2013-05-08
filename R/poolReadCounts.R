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
    newcounts[,x] <- rowSums(counts[,replicates, drop=FALSE])
    if ('corrected' %in% names(obj))
      newcorrected[,x] <- rowSums(corrected[,replicates, drop=FALSE])
    oldpheno <- dat[['phenodata']][replicates,]
    newpheno[x, 'reads'] <- sum(oldpheno$reads)
    if (length(unique(oldpheno$loess.span)) == 1)
      newpheno[x, 'loess.span'] <- oldpheno[1, 'loess.span']
    if (length(unique(oldpheno$loess.family)) == 1)
      newpheno[x, 'loess.family'] <- oldpheno[1, 'loess.family']
  }
  if ('corrected' %in% names(obj))
    return(list(phenodata=newphenodata, bins=bins, counts=newcounts, corrected=newcorrected))
  list(phenodata=newphenodata, bins=bins, counts=newcounts)
}

# EOF

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
  if ('copynumber' %in% assayDataElementNames(object))
    assayDataElement(object, 'copynumber') <- NULL
  if ('segmented' %in% assayDataElementNames(object))
    assayDataElement(object, 'segmented') <- NULL
  if ('calls' %in% assayDataElementNames(object)) {
    assayDataElement(object, 'calls') <- NULL
    assayDataElement(object, 'probloss') <- NULL
    assayDataElement(object, 'probnorm') <- NULL
    assayDataElement(object, 'probgain') <- NULL
    if ('probdloss' %in% assayDataElementNames(object))
      assayDataElement(object, 'probdloss') <- NULL
    if ('probamp' %in% assayDataElementNames(object))
      assayDataElement(object, 'probamp') <- NULL
  }

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
  if (length(newsamples) == nrow(phenodata)) {
    sampleNames(object) <- samples
    object <- object[,order(samples)]
    return(object)
  }
  newphenodata <- data.frame(name=newsamples, reads=NA_integer_, loess.span=NA_real_,
    loess.family=NA_character_, row.names=newsamples, stringsAsFactors=FALSE)
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

setMethod('initialize', 'QDNAseqReadCounts', function(.Object, bins, counts,
  phenodata, ...) {
  callNextMethod(.Object, featureData=AnnotatedDataFrame(bins),
    assayData=assayDataNew(counts=counts),
    phenoData=AnnotatedDataFrame(phenodata), ...)
})

setMethod('binFilter', signature=c(object='QDNAseqReadCounts'),
  definition=function(object) {
  if ('filter' %in% colnames(fData(object)))
    return(fData(object)$filter)
  rep(TRUE, times=nrow(object))
})

setMethod('chromosomes', signature=c(object='QDNAseqReadCounts'),
  definition=function(object) {
  tmp <- fData(object)$chromosome
  tmp[tmp=='X'] <- '23'
  tmp[tmp=='Y'] <- '24'
  as.integer(tmp)
})

setMethod('bpstart', signature=c(object='QDNAseqReadCounts'),
  definition=function(object) {
  fData(object)$start
})

setMethod('bpend', signature=c(object='QDNAseqReadCounts'),
  definition=function(object) {
  fData(object)$end
})

setMethod('copynumber', signature=c(object='QDNAseqReadCounts'),
  definition=function(object) {
  assayDataElement(object, 'copynumber')
})

setMethod('segmented', signature=c(object='QDNAseqReadCounts'),
  definition=function(object) {
  assayDataElement(object, 'segmented')
})

setMethod('calls', signature=c(object='QDNAseqReadCounts'),
  definition=function(object) {
  assayDataElement(object, 'calls')
})

setMethod('probdloss', signature=c(object='QDNAseqReadCounts'),
  definition=function(object) {
  assayDataElement(object, 'probdloss')
})

setMethod('probloss', signature=c(object='QDNAseqReadCounts'),
  definition=function(object) {
  assayDataElement(object, 'probloss')
})

setMethod('probnorm', signature=c(object='QDNAseqReadCounts'),
  definition=function(object) {
  assayDataElement(object, 'probnorm')
})

setMethod('probgain', signature=c(object='QDNAseqReadCounts'),
  definition=function(object) {
  assayDataElement(object, 'probgain')
})

setMethod('probamp', signature=c(object='QDNAseqReadCounts'),
  definition=function(object) {
  assayDataElement(object, 'probamp')
})

setReplaceMethod('binFilter', signature=c(object='QDNAseqReadCounts',
  value='logical'), definition=function(object, value) {
  fData(object)$filter <- value
})

setReplaceMethod('copynumber', signature=c(object='QDNAseqReadCounts',
  value='matrix'), definition=function(object, value) {
  if (nrow(value) == nrow(object)) {
    assayDataElementReplace(object, 'copynumber', value)
  } else {
    value2 <- matrix(nrow=nrow(object), ncol=ncol(object),
      dimnames=list(featureNames(object), sampleNames(object)))
    value2[rownames(value), ] <- value
    assayDataElementReplace(object, 'copynumber', value2)
  }
})

setReplaceMethod('segmented', signature=c(object='QDNAseqReadCounts',
  value='matrix'), definition=function(object, value) {
  if (nrow(value) == nrow(object)) {
    assayDataElementReplace(object, 'segmented', value)
  } else {
    value2 <- matrix(nrow=nrow(object), ncol=ncol(object),
      dimnames=list(featureNames(object), sampleNames(object)))
    value2[rownames(value), ] <- value
    assayDataElementReplace(object, 'segmented', value2)
  }
})

setReplaceMethod('calls', signature=c(object='QDNAseqReadCounts',
  value='matrix'), definition=function(object, value) {
  if (nrow(value) == nrow(object)) {
    assayDataElementReplace(object, 'calls', value)
  } else {
    value2 <- matrix(nrow=nrow(object), ncol=ncol(object),
      dimnames=list(featureNames(object), sampleNames(object)))
    value2[rownames(value), ] <- value
    assayDataElementReplace(object, 'calls', value2)
  }
})

setReplaceMethod('probdloss', signature=c(object='QDNAseqReadCounts',
  value='matrix'), definition=function(object, value) {
  if (nrow(value) == nrow(object)) {
    assayDataElementReplace(object, 'probdloss', value)
  } else {
    value2 <- matrix(nrow=nrow(object), ncol=ncol(object),
      dimnames=list(featureNames(object), sampleNames(object)))
    value2[rownames(value), ] <- value
    assayDataElementReplace(object, 'probdloss', value2)
  }
})

setReplaceMethod('probloss', signature=c(object='QDNAseqReadCounts',
  value='matrix'), definition=function(object, value) {
  if (nrow(value) == nrow(object)) {
    assayDataElementReplace(object, 'probloss', value)
  } else {
    value2 <- matrix(nrow=nrow(object), ncol=ncol(object),
      dimnames=list(featureNames(object), sampleNames(object)))
    value2[rownames(value), ] <- value
    assayDataElementReplace(object, 'probloss', value2)
  }
})

setReplaceMethod('probnorm', signature=c(object='QDNAseqReadCounts',
  value='matrix'), definition=function(object, value) {
  if (nrow(value) == nrow(object)) {
    assayDataElementReplace(object, 'probnorm', value)
  } else {
    value2 <- matrix(nrow=nrow(object), ncol=ncol(object),
      dimnames=list(featureNames(object), sampleNames(object)))
    value2[rownames(value), ] <- value
    assayDataElementReplace(object, 'probnorm', value2)
  }
})

setReplaceMethod('probgain', signature=c(object='QDNAseqReadCounts',
  value='matrix'), definition=function(object, value) {
  if (nrow(value) == nrow(object)) {
    assayDataElementReplace(object, 'probgain', value)
  } else {
    value2 <- matrix(nrow=nrow(object), ncol=ncol(object),
      dimnames=list(featureNames(object), sampleNames(object)))
    value2[rownames(value), ] <- value
    assayDataElementReplace(object, 'probgain', value2)
  }
})

setReplaceMethod('probamp', signature=c(object='QDNAseqReadCounts',
  value='matrix'), definition=function(object, value) {
  if (nrow(value) == nrow(object)) {
    assayDataElementReplace(object, 'probamp', value)
  } else {
    value2 <- matrix(nrow=nrow(object), ncol=ncol(object),
      dimnames=list(featureNames(object), sampleNames(object)))
    value2[rownames(value), ] <- value
    assayDataElementReplace(object, 'probamp', value2)
  }
})

## adapted from CGHcall:

# EOF

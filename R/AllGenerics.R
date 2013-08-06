setGeneric('correctBins', function(object, span=0.65, family='symmetric',
  adjustIncompletes=TRUE, keepCounts=TRUE, storeResiduals=TRUE, force=FALSE,
  ...)
  standardGeneric('correctBins'))
setGeneric('poolRuns', function(object, samples, force=FALSE)
  standardGeneric('poolRuns'))
setGeneric('compareToReference', function(object, references, force=FALSE)
  standardGeneric('compareToReference'))
setGeneric('normalizeBins', function(object, method='median',
  smoothOutliers=TRUE, logTransform=TRUE, force=FALSE, ...)
  standardGeneric('normalizeBins'))
setGeneric('highlightFilters', function(object, col='red', mappability=50,
  blacklist=0, residual=1, bases=100, ...)
  standardGeneric('highlightFilters'))
setGeneric('applyFilters', function(object, mappability=50,
  blacklist=0, residual=1, bases=100, filterAllosomes=TRUE, force=FALSE)
  standardGeneric('applyFilters'))
setGeneric('segmentBins', function(object, weights=FALSE, normalize=TRUE,
  inter=c(-0.1,0.1), force=FALSE, ...) standardGeneric('segmentBins'))
setGeneric('callBins', function(object, ...)
  standardGeneric('callBins'))
setGeneric('makeCgh', function(object) standardGeneric('makeCgh'))

setGeneric('readCountPlot', function(x, y, ...)
  standardGeneric('readCountPlot'))






setGeneric('binFilter', function(object) standardGeneric('binFilter'))

setGeneric('chromosomes', function(object) standardGeneric('chromosomes'))
setGeneric('bpstart', function(object) standardGeneric('bpstart'))
setGeneric('bpend', function(object) standardGeneric('bpend'))

setGeneric('copynumber', function(object) standardGeneric('copynumber'))
setGeneric('segmented', function(object) standardGeneric('segmented'))
setGeneric('calls', function(object) standardGeneric('calls'))
setGeneric('probdloss', function(object) standardGeneric('probdloss'))
setGeneric('probloss', function(object) standardGeneric('probloss'))
setGeneric('probnorm', function(object) standardGeneric('probnorm'))
setGeneric('probgain', function(object) standardGeneric('probgain'))
setGeneric('probamp', function(object) standardGeneric('probamp'))

setGeneric('binFilter<-', function(object, value)
  standardGeneric('binFilter<-'))

setGeneric('copynumber<-', function(object, value)
  standardGeneric('copynumber<-'))
setGeneric('segmented<-', function(object, value)
  standardGeneric('segmented<-'))
setGeneric('calls<-', function(object, value)
  standardGeneric('calls<-'))
setGeneric('probdloss<-', function(object, value)
  standardGeneric('probdloss<-'))
setGeneric('probloss<-', function(object, value)
  standardGeneric('probloss<-'))
setGeneric('probnorm<-', function(object, value)
  standardGeneric('probnorm<-'))
setGeneric('probgain<-', function(object, value)
  standardGeneric('probgain<-'))
setGeneric('probamp<-', function(object, value)
  standardGeneric('probamp<-'))

# EOF
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
setGeneric('highlightFilters', function(object, col='red', residual=NA,
  blacklist=NA, mappability=NA, bases=NA, type='union', ...)
  standardGeneric('highlightFilters'))
setGeneric('applyFilters', function(object, residual=4, blacklist=0,
  mappability=50, bases=99, filterAllosomes=TRUE, force=FALSE)
  standardGeneric('applyFilters'))
setGeneric('segmentBins', function(object, weights=FALSE, normalize=TRUE,
  inter=c(-0.1,0.1), force=FALSE, ...) standardGeneric('segmentBins'))
setGeneric('callBins', function(object, ...)
  standardGeneric('callBins'))
setGeneric('makeCgh', function(object) standardGeneric('makeCgh'))

setGeneric('readCountPlot', function(x, y, ...)
  standardGeneric('readCountPlot'))
setGeneric('noisePlot', function(x, y, ...)
  standardGeneric('noisePlot'))






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
#setGeneric('plot', function(x, y) standardGeneric('plot'))

setGeneric('applyFilters', function(object, mappability=50, blacklist=0,
  residual=1, bases=100, filterAllosomes=TRUE, force=FALSE)
   standardGeneric('applyFilters'))
setGeneric('CGHcall', function(object, ...)
  standardGeneric('CGHcall'))
setGeneric('compareToReference', function(object, references)
  standardGeneric('compareToReference'))
setGeneric('correct', function(object, span=0.65, family='symmetric',
  adjustIncompletes=TRUE, keepCounts=TRUE, storeResiduals=TRUE, ...)
  standardGeneric('correct'))
setGeneric('ExpandCGHcall', function(listcall, object, ...)
  standardGeneric('ExpandCGHcall'))
setGeneric('highlightFilters', function(object, col='red', mappability=50, blacklist=0, residual=1, bases=100, ...) standardGeneric('highlightFilters'))
setGeneric('makeCgh', function(object) standardGeneric('makeCgh'))
setGeneric('normalize', function(object, method='median', smoothOutliers=TRUE,
  logTransform=TRUE, ...) standardGeneric('normalize'))
setGeneric('poolRuns', function(object, samples) standardGeneric('poolRuns'))
setGeneric('postsegnormalize', function(object, ...)
  standardGeneric('postsegnormalize'))
setGeneric('segmentData', function(object, weights=TRUE, ...)
  standardGeneric('segmentData'))







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
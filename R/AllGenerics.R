setGeneric("applyFilters", function(object, residual=TRUE, blacklist=TRUE,
  mappability=NA, bases=NA, filterAllosomes=TRUE, force=FALSE)
  standardGeneric("applyFilters"))
setGeneric("callBins", function(object, ...)
  standardGeneric("callBins"))
setGeneric("compareToReference", function(object, references, offset=2^-10,
  force=FALSE)
  standardGeneric("compareToReference"))
setGeneric("correctBins", function(object, fit=NULL,
  type="ratio", adjustIncompletes=TRUE, ...)
  standardGeneric("correctBins"))
setGeneric("estimateCorrection", function(object, span=0.65, family="symmetric",
  adjustIncompletes=TRUE, maxIter=1, cutoff=4.0, ...)
  standardGeneric("estimateCorrection"))
setGeneric("highlightFilters", function(object, col="red", residual=NA,
  blacklist=NA, mappability=NA, bases=NA, type="union",
  logTransform=TRUE, ...)
  standardGeneric("highlightFilters"))
setGeneric("isobarPlot", function(x, y, ...)
  standardGeneric("isobarPlot"))
setGeneric("makeCgh", function(object, filter=TRUE, ...)
  standardGeneric("makeCgh"))
setGeneric("noisePlot", function(x, y, ...)
  standardGeneric("noisePlot"))
setGeneric("normalizeBins", function(object, method="median", force=FALSE)
  standardGeneric("normalizeBins"))
setGeneric("normalizeSegmentedBins", function(object, inter=c(-0.1, 0.1),
  force=FALSE, ...) standardGeneric("normalizeSegmentedBins"))
setGeneric("poolRuns", function(object, samples, force=FALSE)
  standardGeneric("poolRuns"))
setGeneric("segmentBins", function(object, smoothBy=FALSE,
  alpha=1e-10, undo.splits="sdundo", undo.SD=1.0,
  force=FALSE, ...) standardGeneric("segmentBins"))
setGeneric("smoothOutlierBins", function(object,
  logTransform=TRUE, force=FALSE, ...)
  standardGeneric("smoothOutlierBins"))


setGeneric("binsToUse", function(object) standardGeneric("binsToUse"))

setGeneric("chromosomes", function(object) standardGeneric("chromosomes"))
setGeneric("bpstart", function(object) standardGeneric("bpstart"))
setGeneric("bpend", function(object) standardGeneric("bpend"))

setGeneric("counts", function(object) standardGeneric("counts"))
setGeneric("fit", function(object) standardGeneric("fit"))
setGeneric("copynumber", function(object) standardGeneric("copynumber"))
setGeneric("segmented", function(object) standardGeneric("segmented"))
setGeneric("calls", function(object) standardGeneric("calls"))
setGeneric("probdloss", function(object) standardGeneric("probdloss"))
setGeneric("probloss", function(object) standardGeneric("probloss"))
setGeneric("probnorm", function(object) standardGeneric("probnorm"))
setGeneric("probgain", function(object) standardGeneric("probgain"))
setGeneric("probamp", function(object) standardGeneric("probamp"))

setGeneric("binsToUse<-", function(object, value)
  standardGeneric("binsToUse<-"))

setGeneric("counts<-", function(object, value)
  standardGeneric("counts<-"))
setGeneric("fit<-", function(object, value)
  standardGeneric("fit<-"))
setGeneric("copynumber<-", function(object, value)
  standardGeneric("copynumber<-"))
setGeneric("segmented<-", function(object, value)
  standardGeneric("segmented<-"))
setGeneric("calls<-", function(object, value)
  standardGeneric("calls<-"))
setGeneric("probdloss<-", function(object, value)
  standardGeneric("probdloss<-"))
setGeneric("probloss<-", function(object, value)
  standardGeneric("probloss<-"))
setGeneric("probnorm<-", function(object, value)
  standardGeneric("probnorm<-"))
setGeneric("probgain<-", function(object, value)
  standardGeneric("probgain<-"))
setGeneric("probamp<-", function(object, value)
  standardGeneric("probamp<-"))

# EOF

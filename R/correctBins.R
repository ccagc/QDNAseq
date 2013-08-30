#########################################################################/**
# @RdocFunction correctBins
#
# @alias correctBins,QDNAseqReadCounts-method
#
# @title "Correct binned read counts for GC content and mappability"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{object}{...}
#   \item{span}{...}
#   \item{family}{...}
#   \item{adjustIncompletes}{...}
#   \item{keepCounts}{...}
#   \item{storeResiduals}{...}
#   \item{...}{Further aguments to loess.}
#   \item{force}{...}
# }
#
# \value{
#   Returns a named @list with elements ...
# }
#
# @author "IS"
#
# \seealso{
#   Internally, @see "stats::loess" is used to fit the regression model.
# }
#*/#########################################################################

setMethod('correctBins', signature=c(object='QDNAseqReadCounts'),
  definition=function(object, span=0.65, family='symmetric',
  adjustIncompletes=TRUE, keepCounts=TRUE, storeResiduals=FALSE, force=FALSE,
  ...) {
  if (!force && 'copynumber' %in% assayDataElementNames(object))
    stop('Data has already been normalized. Changing the correction will ',
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
  counts <- assayDataElement(object, 'counts')
  if (adjustIncompletes) {
    counts <- counts / fData(object)$bases * 100L
    counts[fData(object)$bases == 0] <- 0L
  }
  message('Performing correction for GC content and mappability:')
  if (length(span) == 1L)
    span <- rep(span, times=ncol(counts))
  if (length(family) == 1L)
    family <- rep(family, times=ncol(counts))
  if (length(span) != ncol(counts))
    stop('Parameter span has to be either a single value or a vector the ',
      'same length as there are samples in object.')
  if (length(family) != ncol(counts))
    stop('Parameter family has to be either a single value or a vector the ',
      'same length as there are samples in object.')
  if ('filter' %in% colnames(fData(object))) {
    condition <- fData(object)$filter
  } else {
    condition <- rep(TRUE, times=nrow(object))
  }
  used.span <- rep(NA_real_, times=ncol(counts))
  used.family <- rep(NA_real_, times=ncol(counts))
  corrected <- matrix(nrow=nrow(counts), ncol=ncol(counts),
    dimnames=dimnames(counts))
  residuals <- matrix(nrow=nrow(counts), ncol=ncol(counts),
    dimnames=dimnames(counts))
  gc <- round(fData(object)$gc)
  mappability <- round(fData(object)$mappability)
  median.counts <- aggregate(counts[condition, ], by=list(gc=gc[condition],
    mappability=mappability[condition]), FUN=median)
  median.counts <- median.counts[!is.na(median.counts$gc), ]
  median.counts <- median.counts[!is.na(median.counts$mappability), ]
  rownames(median.counts) <- paste(median.counts$gc, '-',
    median.counts$mappability, sep='')
  all.combinations <- expand.grid(gc=unique(gc[!is.na(gc)]),
    mappability=unique(mappability[!is.na(mappability)]))
  rownames(all.combinations) <- paste(all.combinations$gc, '-',
    all.combinations$mappability, sep='')
  for (i in seq_len(ncol(counts))) {
    if (is.na(span[i]) && is.na(family[i])) {
      message('  Skipping correction for sample ', colnames(counts)[i], '...')
      next
    }
    message('  Using span=', span[i], '\tand family=', family[i],
      ',\tcorrecting sample ', colnames(counts)[i], '...')
    vals <- median.counts[, i+2L]
    corvals <- counts[, i]
    try({
      l <- loess(vals ~ gc * mappability, data=median.counts,
        span=span[i], family=family[i]) # , ...)
      fit <- as.vector(predict(l, all.combinations))
      names(fit) <- rownames(all.combinations)
      residuals[, i] <- (corvals - fit[paste(gc, '-', mappability, sep='')]) /
        fit[paste(gc, '-', mappability, sep='')]
      correction <- median(fit, na.rm=TRUE) - fit
      corvals <- corvals + correction[paste(gc, '-', mappability, sep='')]
      corvals <- corvals - min(corvals, na.rm=TRUE)
      used.span[i] <- span[i]
      used.family[i] <- family[i]
    }, silent=TRUE)
    corrected[, i] <- corvals
  }
  corrected <- round(corrected, digits=2L)
  object$loess.span <- used.span
  object$loess.family <- used.family
  assayDataElement(object, 'corrected') <- corrected
  if (!keepCounts)
    assayDataElement(object, 'counts') <- NULL
  if (storeResiduals)
    assayDataElement(object, 'residuals') <- residuals
  message('Done.')
  object
})

# a subfunction for performing the correction on one sample at a time could
# be defined here.

# EOF

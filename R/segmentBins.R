#########################################################################/**
# @RdocFunction segmentBins
#
# @alias segmentBins,QDNAseqCopyNumbers-method
#
# @title "Segments normalized copy number data"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{object}{An object of class QDNAseqCopyNumbers.}
#   \item{smoothBy}{An optional integer value to perform smoothing before
#     segmentation by taking the mean of every smoothBy bins, and then segment
#     those means. Default is to perform no smoothing.}
#   \item{alpha}{Significance levels for the test to accept change-points.}
#   \item{undo.splits}{A character string specifying how change-points are to be
#     undone, if at all.  Default is "none".  Other choices are
#     "prune", which uses a sum of squares criterion, and "sdundo",
#     which undoes splits that are not at least this many SDs
#     apart.}
#   \item{undo.SD}{The number of SDs between means to keep a split if
#     undo.splits="sdundo".}
#   \item{normalize}{Whether to perform normalization after segmentation with
#     @see "CGHcall::postsegnormalize" of the \pkg{CGHcall} package.}
#   \item{inter}{If performing normalization, the interval in which the
#     function should search for the normal level.}
#   \item{force}{Whether to force execution when it causes removal of
#     downstream calling results.}
#   \item{...}{Additional arguments passed to @see "DNAcopy::segment".}
# }
#
# \value{
#   Returns an object of class QDNAseqCopyNumbers with segmentation results
#     added.
# }
#
# @author "IS"
#
# \seealso{
#   Internally, @see "DNAcopy::segment" of the \pkg{DNAcopy} package,
#   which implements the CBS method, is used to segment the data.
#   Optionally, normalization is performed with @see "CGHcall::postsegnormalize"
#   of the \pkg{CGHcall} package.
# }
#
#*/#########################################################################
setMethod('segmentBins', signature=c(object='QDNAseqCopyNumbers'),
  definition=function(object, smoothBy=FALSE, alpha=1e-10,
  undo.splits='sdundo', undo.SD=1.0,
  normalize=TRUE, inter=c(-0.1, 0.1), force=FALSE, ...) {

  if (!force && 'calls' %in% assayDataElementNames(object))
    stop('Data has already been called. Changing the segmentation will ',
      'remove calling ',
      'results. Please specify force=TRUE, if you want this.')
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
  condition <- binsToUse(object)

  if (!is.numeric(smoothBy) || smoothBy <= 1) {
    vmsg('Performing segmentation:')
  } else {
    vmsg('Performing segmentation with smoothing over ',
      smoothBy, ' bins:')
  }

  copynumber <- copynumber(object)
  copynumber[!condition, ] <- NA_real_
  copynumber <- log2adhoc(copynumber)
  segmented <- matrix(NA_real_, nrow=nrow(copynumber), ncol=ncol(copynumber),
    dimnames=dimnames(copynumber))

  ## loop through samples
  for (s in seq_len(ncol(copynumber))) {
    vmsg('  Segmenting: ', sampleNames(object)[s],
      ' (', s, ' of ', ncol(object), ') ...', appendLF=FALSE)

    ## loop through chromosomes
    for (chr in unique(fData(object)$chromosome[condition])) {
      index <- fData(object)$chromosome == chr
      chrStarts <- fData(object)$start[index]
      ## smooth if needed
      if (is.na(smoothBy) || !smoothBy || smoothBy <= 1) {
        chrCopynumber <- copynumber[index, s]
      } else {
        binToBin <- 0:(sum(index)-1) %/% smoothBy
        chrCopynumber <- aggregate(copynumber[index, s],
          by=list(binToBin), mean, na.rm=TRUE)$x
        chrStarts <- chrStarts[seq(from=1, by=smoothBy,
          length.out=length(chrCopynumber))]
      }

      ## segment
      cna <- CNA(genomdat=chrCopynumber,
        chrom=chr, maploc=chrStarts, data.type='logratio',
        presorted=TRUE)
      segments <- segment(cna, verbose=0,
        alpha=alpha, undo.splits=undo.splits, undo.SD=undo.SD, ...)
      chrSegmented <- rep(NA_real_, length=length(chrCopynumber))
      for (i in 1:nrow(segments$output))
        chrSegmented[segments$segRows$startRow[i]:
          segments$segRows$endRow[i]] <- segments$output$seg.mean[i]

      ## process results whether smoothed or not
      if (is.na(smoothBy) || !smoothBy || smoothBy <= 1) {
        segmented[index, s] <- chrSegmented
      } else {
        segmented[index, s] <- rep(chrSegmented, times=table(binToBin))
      }
    }
    vmsg()
  }
  segmented[is.na(copynumber)] <- NA_real_

  if (!normalize) {
    segmented <- unlog2adhoc(segmented)
    segmented(object) <- segmented
    return(object)
  }

  ## adapted from CGHcall::postsegnormalize()
  seg <- segmented
  values <- colMedians(seg, na.rm=TRUE)
  seg <- t(t(seg) - values)
  countlevall <- apply(seg, MARGIN=2L, FUN=function(x)
    as.data.frame(table(x)))

  intcount <- function(int, sv){
    sv1 <- as.numeric(as.vector(sv[, 1L]))
    wh <- which(sv1 <= int[2L] & sv1 >= int[1L])
    return(sum(sv[wh, 2L]))
  }

  postsegnorm <- function(segvec, int=inter, intnr=3){
    intlength <- (int[2L]-int[1L])/2
    gri <- intlength/intnr
    intst <- int[1L]+(0:intnr)*gri
    intend <- intst+intlength
    ints <- cbind(intst, intend)
    intct <- apply(ints, MARGIN=1L, FUN=intcount, sv=segvec)
    whmax <- which.max(intct)
    return(ints[whmax, ])
  }

  postsegnorm_rec <- function(segvec, int, intnr=3){
    newint <- postsegnorm(segvec, int, intnr)
    newint <- postsegnorm(segvec, newint, intnr)
    newint <- postsegnorm(segvec, newint, intnr)
    newint <- postsegnorm(segvec, newint, intnr)
    newint <- postsegnorm(segvec, newint, intnr)
    return(newint[1L]+(newint[2L]-newint[1L])/2)
  }

  listres <- lapply(countlevall, FUN=postsegnorm_rec, int=inter)
  vecres <- c()
  for (i in seq_along(listres))
    vecres <- c(vecres, listres[[i]])

  segmented <- t(t(seg) - vecres)
  segmented <- unlog2adhoc(segmented)
  segmented[segmented < 0] <- 0
  segmented(object) <- segmented
  copynumber <- t(t(copynumber) - values - vecres)
  copynumber <- unlog2adhoc(copynumber)
  copynumber[copynumber < 0] <- 0
  copynumber(object) <- copynumber
  object
})

# EOF

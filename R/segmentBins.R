#########################################################################/**
# @RdocFunction segmentBins
#
# @alias segmentBins,QDNAseqReadCounts-method
#
# @title "Segments and calls total copy numbers"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{object}{...}
#   \item{smoothBy}{...}
#   \item{alpha}{...}
#   \item{undo.splits}{...}
#   \item{undo.SD}{...}
#   \item{normalize}{...}
#   \item{inter}{...}
#   \item{force}{...}
#   \item{...}{Additional arguments passed to @see "DNAcopy::segment".}
# }
#
# \value{
#   Returns ...
# }
#
# @author "IS"
#
# \seealso{
#   Internally, @see "DNAcopy::segment" of the \pkg{DNAcopy} package,
#   which implements the CBS method, is used to segment the data.
# }
#
#*/#########################################################################
setMethod('segmentBins', signature=c(object='QDNAseqReadCounts'),
  definition=function(object, smoothBy=FALSE, alpha=1e-10,
  undo.splits='sdundo', undo.SD=1.0, normalize=TRUE,
  inter=c(-0.1, 0.1), force=FALSE, ...) {

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
  if ('filter' %in% colnames(fData(object))) {
    condition <- fData(object)$filter
  } else {
    condition <- rep(TRUE, times=nrow(object))
  }

  if (is.na(smoothBy) || !smoothBy || smoothBy <= 1) {
    message('Performing segmentation:')
  } else {
    message('Performing segmentation with smoothing over ',
      smoothBy, ' bins:')
  }

  copynumber <- copynumber(object)
  copynumber[!condition, ] <- NA_real_
  segmented <- matrix(NA_real_, nrow=nrow(copynumber), ncol=ncol(copynumber),
    dimnames=dimnames(copynumber))

  ## loop through samples
  for (s in seq_len(ncol(copynumber))) {
    message('  Segmenting: ', sampleNames(object)[s],
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
  message()
  }
  segmented[is.na(copynumber)] <- NA_real_

  if (!normalize) {
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

  segmented(object) <- t(t(seg) - vecres)
  copynumber(object) <- t(t(copynumber) - values - vecres)
  object
})

# EOF

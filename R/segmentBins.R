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
#   \item{weights}{Either @TRUE or a vector of weights. If @TRUE,
#     loess residuals are used as weights.}
#   \item{normalize}{...}
#   \item{inter}{...}
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
  definition=function(object, weights=TRUE, normalize=TRUE,
  inter=c(-0.1,0.1), ...) {
  if (length(weights) == 1L & weights) {
    if ('residual' %in% colnames(fData(object))) {
      message('Using median loess residuals of control data set as ',
        'segmentation weights.')
      residual <- fData(object)$residual
    } else if ('residuals' %in% assayDataElementNames(object)) {
      message('Using median loess residuals as segmentation weights.')
      residual <- rowMedians(assayDataElement(object, 'residuals'), na.rm=TRUE)
    } else {
      stop('No loess residuals found. Please provide a vector of weights ',
        'or specify weights=FALSE.')
    }
    if (any(is.na(residual[binFilter(object)]))) {
      message('Filtering out ', sum(is.na(residual[binFilter(object)])),
        ' bins with missing residuals.')
      binFilter(object) <- binFilter(object) & !is.na(residual)
    }
    residual <- residual[binFilter(object)]
    residual <- abs(residual)
    residual[residual == 0] <- min(residual[residual != 0], na.rm=TRUE)
    weights <- 1/residual
  }
  copynumber <- copynumber(object)[binFilter(object), , drop=FALSE]
  CNA.object <- CNA(genomdat=copynumber,
    chrom=chromosomes(object)[binFilter(object)],
    maploc=bpstart(object)[binFilter(object)], data.type="logratio",
    sampleid=paste(sampleNames(object), ':', 1:ncol(object), 'of',
      ncol(object), sep=''), presorted=TRUE)
  message('Start data segmentation ...')
  if (length(weights) == 1L && !weights) {
    segmented <- segment(CNA.object, ...)
  } else {
    segmented <- segment(CNA.object, weights=weights, ...)
  }
  numclone <- segmented$output$num.mark
  smrat <- segmented$output$seg
  numsmrat <- cbind(smrat, numclone)
  repdata <- function(row) {
    rep(row[1L], times=row[2L])
  }
  makelist <- apply(numsmrat, MARGIN=1L, FUN=repdata)
  joined <- unlist(makelist)
  ## Not needed anymore
  rm(list="makelist")
  joined <- matrix(joined, ncol=ncol(object), byrow=FALSE)
  rownames(joined) <- rownames(copynumber)
  colnames(joined) <- colnames(copynumber)

  if (!normalize) {
    segmented(object) <- joined
    return(object)
  }

  ## adapted from CGHcall::postsegnormalize()
  seg <- joined
  values <- colMedians(seg, na.rm=TRUE)
  seg <- t(t(seg) - values)
  countlevall <- apply(seg, MARGIN=2L, FUN=function(x) as.data.frame(table(x)))

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

  listres <- lapply(countlevall, postsegnorm_rec, int=inter)
  vecres <- c()
  for(i in 1:length(listres))
    vecres <- c(vecres, listres[[i]])

  segmented(object) <- t(t(seg) - vecres)
  copynumber(object) <- t(t(copynumber) - values - vecres)
  object
})

# EOF

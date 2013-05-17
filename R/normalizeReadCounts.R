normalizeReadCounts <- function(obj, method='median', logTransform=TRUE, smoothOutliers=TRUE) {
  if (exists('filter', obj)) {
    condition <- obj[['filter']]
  } else {
    condition <- rep(TRUE, nrow(obj[['bins']]))
  }
  bins <- obj[['bins']][condition,]
  copynumber <- obj[['corrected']][condition,]
  if (logTransform)
    copynumber <- log2(copynumber + 1)
  if (method == 'none') {
    cat('Skipping normalization ... \n')
  } else {
    if (method == 'median') {
      cat('Applying median normalization ... \n')
      values <- apply(copynumber, 2, median, na.rm=TRUE)
    } else if (method == 'mode') {
      cat('Applying mode normalization ... \n')
      values <- apply(copynumber, 2, function(x) { d <- density(x, na.rm=TRUE); d$x[which.max(d$y)] })
    }
    copynumber <- t(t(copynumber) - values)
  }
  if (smoothOutliers) {
    cat('Smoothing outliers ... \n')
    CNA.object <- DNAcopy::smooth.CNA(DNAcopy::CNA(copynumber, bins$chromosome, bins$start, data.type='logratio', presorted=TRUE))
    copynumber <- as.matrix(CNA.object[,-(1:2)])
  }
  obj[['copynumber']] <- matrix(nrow=nrow(obj[['corrected']]), ncol=ncol(obj[['corrected']]), dimnames=dimnames(obj[['corrected']]))
  obj[['copynumber']][rownames(copynumber),] <- copynumber
  obj
}

# EOF

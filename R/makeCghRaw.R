makeCghRaw <- function(obj) {
  if (exists('filter', obj)) {
    condition <- obj[['filter']]
  } else {
    condition <- rep(TRUE, nrow(obj[['bins']]))
  }
  cgh <- make_cghRaw(data.frame(bin=rownames(obj[['bins']][condition,]), obj[['bins']][condition, c('chromosome', 'start', 'end'),], obj[['copynumber']][condition,], check.names=FALSE, stringsAsFactors=FALSE))
  pData(cgh) <- obj[['phenodata']]
  cgh
}

# EOF

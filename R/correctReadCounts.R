correctReadCounts <- function(dat, span=0.65, family='symmetric', plotting=FALSE, blacklist.threshold=0, ...) {
  if (length(span) != 1 && length(span) != ncol(dat[['counts']]))
    stop('span has to be either a single value or a vector the same length as the number of samples in dat.')
  if (length(family) != 1 && length(family) != ncol(dat[['counts']]))
    stop('family has to be either a single value or a vector the same length as the number of samples in dat.')
  if (length(span) == 1)
    span <- rep(span, ncol(dat[['counts']]))
  if (length(family) == 1)
    family <- rep(family, ncol(dat[['counts']]))
  used.span <- rep(NA, ncol(dat[['counts']]))
  used.family <- rep(NA, ncol(dat[['counts']]))
  cat('Performing correction for GC content and mappability:\n')
  dat[['corrected']] <- matrix(nrow=nrow(dat[['counts']]), ncol=ncol(dat[['counts']]), dimnames=dimnames(dat[['counts']]))
  dat[['residuals']] <- matrix(nrow=nrow(dat[['counts']]), ncol=ncol(dat[['counts']]), dimnames=dimnames(dat[['counts']]))
  condition <- dat[['bins']]$chromosome %in% as.character(1:22) & dat[['bins']]$blacklist <= blacklist.threshold
  gc <- round(dat[['bins']]$gc)
  mappability <- round(dat[['bins']]$mappability)
  median.counts <- aggregate(dat[['counts']][condition,], by=list(gc=gc[condition], mappability=mappability[condition]), median)
  median.counts <- median.counts[!is.na(median.counts$gc),]
  rownames(median.counts) <- paste(median.counts$gc, '-', median.counts$mappability, sep='')
  residuals <- matrix(nrow=nrow(median.counts), ncol=ncol(dat[['counts']]), dimnames=list(rownames(median.counts), colnames(dat[['counts']])))
  if (plotting) {
    x <- min(median.counts$mappability):max(median.counts$mappability)
    y <- min(median.counts$gc):max(median.counts$gc)
    m <- matrix(nrow=length(x), ncol=length(y), dimnames=list(x, y))
    d <- sub('-counts\\.(txt|RData)', '-correction', basename(countsfile))
    if (!file.exists(d))
      dir.create(d)
    setwd(d)
  }
  for (i in 1:ncol(dat[['counts']])) {
    if (is.na(span[i]) && is.na(family[i])) {
      cat('\tSkipping correction for sample ', colnames(dat[['counts']])[i], '...\n', sep='')
      next
    }
    cat('\tUsing span=', span[i], '\tand family=', family[i], ',\tcorrecting sample ', colnames(dat[['counts']])[i], '...\n', sep='')
    vals <- median.counts[,i+2]
    corvals <- dat[['counts']][,i]
    try({
      l <- loess(vals ~ median.counts$gc * median.counts$mappability, span=span[i], family=family[i]) # , ...)
      fit <- l$fitted
      names(fit) <- rownames(median.counts)
      dat[['residuals']][,i] <- corvals - fit[paste(gc, '-', mappability, sep='')]
      correction <- median(fit, na.rm=TRUE) - fit
      corvals <- corvals + correction[paste(gc, '-', mappability, sep='')]
      corvals <- corvals - min(corvals, na.rm=TRUE)
      used.span[i] <- span[i]
      used.family[i] <- family[i]
      residuals[,i] <- l$residuals
      if (plotting) {
        pnga4(paste(colnames(dat[['counts']])[i], '-counts.png', sep=''))
        for (j in 1:nrow(median.counts))
          m[as.character(median.counts[j, 'mappability']), as.character(median.counts[j, 'gc'])] <- median.counts[j, i+2]
        image(x, y, m, col=paste('#', c(sprintf('%02X', 0:255), rep('FF', 256)), c(rep('FF', 256), sprintf('%02X', 255:0)), sprintf('%02X', 255), sep=''), xlab='mappability', ylab='GC content', main=paste('Median Read Counts -', colnames(dat[['counts']])[i]), zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))))
        contour(x, y, m, nlevels=20, zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))), add=TRUE)
        dev.off()
        pnga4(paste(colnames(dat[['counts']])[i], '-fit.png', sep=''))
        for (j in 1:nrow(median.counts))
          m[as.character(median.counts[j, 'mappability']), as.character(median.counts[j, 'gc'])] <- fit[j]
        image(x, y, m, col=paste('#', c(sprintf('%02X', 0:255), rep('FF', 256)), c(rep('FF', 256), sprintf('%02X', 255:0)), sprintf('%02X', 255), sep=''), xlab='mappability', ylab='GC content', main=paste('Loess Fit -', colnames(dat[['counts']])[i]), zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))))
        contour(x, y, m, nlevels=20, zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))), add=TRUE)
        dev.off()
        pnga4(paste(colnames(dat[['counts']])[i], '-residuals.png', sep=''))
        for (j in 1:nrow(median.counts))
          m[as.character(median.counts[j, 'mappability']), as.character(median.counts[j, 'gc'])] <- l$residuals[j]
        image(x, y, m, col=paste('#', c(sprintf('%02X', 0:255), rep('FF', 256)), c(rep('FF', 256), sprintf('%02X', 255:0)), sprintf('%02X', 255), sep=''), xlab='mappability', ylab='GC content', main=paste('Loess Residuals -', colnames(dat[['counts']])[i]), zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))))
        contour(x, y, m, nlevels=20, zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))), add=TRUE)
        dev.off()
      }
    }, silent=TRUE)
    dat[['corrected']][,i] <- corvals
  }
  dat[['corrected']] <- round(dat[['corrected']], digits=2)
  median.residuals <- round(apply(residuals, 1, median, na.rm=TRUE), digits=2)
  # dat[['bins']]$median.residual <- median.residuals[paste(gc, '-', mappability, sep='')]
  # dat[['bins']]$median.residual <- round(apply(dat[['residuals']], 1, median, na.rm=TRUE), digits=2)
  if (plotting) {
    for (j in 1:nrow(median.counts))
      m[as.character(median.counts[j, 'mappability']), as.character(median.counts[j, 'gc'])] <- median.residuals[j]
    setwd('..')
    pnga4(sub('-counts\\.(txt|RData)', '-residuals.png', basename(countsfile)))
    image(x, y, m, col=paste('#', c(sprintf('%02X', 0:255), rep('FF', 256)), c(rep('FF', 256), sprintf('%02X', 255:0)), sprintf('%02X', 255), sep=''), xlab='mappability', ylab='GC content', main='Median Residuals', zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))))
    contour(x, y, m, nlevels=20, zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))), add=TRUE)
    dev.off()
  }
  dat[['phenodata']]$loess.span <- used.span
  dat[['phenodata']]$loess.family <- used.family
  cat('Done.\n')
  dat
}

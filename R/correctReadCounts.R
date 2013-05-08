correctReadCounts <- function(obj, span=0.65, family='symmetric', plotting=FALSE, condition=TRUE, ...) {
  phenodata <- obj[['phenodata']]
  bins <- obj[['bins']]
  counts <- obj[['counts']]
  cat('Performing correction for GC content and mappability:\n')
  if (length(span) == 1)
    span <- rep(span, ncol(counts))
  if (length(family) == 1)
    family <- rep(family, ncol(counts))
  if (length(span) != ncol(counts))
    stop('span has to be either a single value or a vector the same length as the number of columns in counts.')
  if (length(family) != ncol(counts))
    stop('family has to be either a single value or a vector the same length as the number of columns in counts.')
  if (length(condition)==1) {
    if (condition) {                                         # when calculating the correction, ignore
      condition <- bins$chromosome %in% as.character(1:22) & # sex chromosomes
        bins$chromosome == c(bins$chromosome[-1], '0')       # and the very last bins of each chromosome,
      if ('blacklist' %in% colnames(bins))                   # and if the information is present,
        condition <- condition & bins$blacklist == 0         # bins overlapping with the ENCDOE blacklist
    } else {
      condition <- rep(TRUE, nrow(bins))
    }
  }
  if (length(condition) != nrow(bins))
    stop('condition has to be either a single value or vector the same length as the number of rows in bins.')
  used.span <- rep(NA, ncol(counts))
  used.family <- rep(NA, ncol(counts))
  corrected <- matrix(nrow=nrow(counts), ncol=ncol(counts), dimnames=dimnames(counts))
  residuals <- matrix(nrow=nrow(counts), ncol=ncol(counts), dimnames=dimnames(counts))
  gc <- round(bins$gc)
  mappability <- round(bins$mappability)
  median.counts <- aggregate(counts[condition], by=list(gc=gc[condition], mappability=mappability[condition]), median)
  median.counts <- median.counts[!is.na(median.counts$gc),]
  rownames(median.counts) <- paste(median.counts$gc, '-', median.counts$mappability, sep='')
  # if (plotting) {
    # x <- min(median.counts$mappability):max(median.counts$mappability)
    # y <- min(median.counts$gc):max(median.counts$gc)
    # m <- matrix(nrow=length(x), ncol=length(y), dimnames=list(x, y))
    # d <- sub('-counts\\.(txt|RData)', '-correction', basename(countsfile))
    # if (!file.exists(d))
      # dir.create(d)
    # setwd(d)
  # }
  for (i in 1:ncol(counts)) {
    if (is.na(span[i]) && is.na(family[i])) {
      cat('\tSkipping correction for sample ', colnames(counts)[i], '...\n', sep='')
      next
    }
    cat('\tUsing span=', span[i], '\tand family=', family[i], ',\tcorrecting sample ', colnames(counts)[i], '...\n', sep='')
    vals <- median.counts[,i+2]
    corvals <- counts[,i]
    try({
      l <- loess(vals ~ median.counts$gc * median.counts$mappability, span=span[i], family=family[i]) # , ...)
      fit <- l$fitted
      names(fit) <- rownames(median.counts)
      residuals[,i] <- corvals - fit[paste(gc, '-', mappability, sep='')]
      correction <- median(fit, na.rm=TRUE) - fit
      corvals <- corvals + correction[paste(gc, '-', mappability, sep='')]
      corvals <- corvals - min(corvals, na.rm=TRUE)
      used.span[i] <- span[i]
      used.family[i] <- family[i]
      # if (plotting) {
        # pnga4(paste(colnames(counts)[i], '-counts.png', sep=''))
        # for (j in 1:nrow(median.counts))
          # m[as.character(median.counts[j, 'mappability']), as.character(median.counts[j, 'gc'])] <- median.counts[j, i+2]
        # image(x, y, m, col=paste('#', c(sprintf('%02X', 0:255), rep('FF', 256)), c(rep('FF', 256), sprintf('%02X', 255:0)), sprintf('%02X', 255), sep=''), xlab='mappability', ylab='GC content', main=paste('Median Read Counts -', colnames(counts)[i]), zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))))
        # contour(x, y, m, nlevels=20, zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))), add=TRUE)
        # dev.off()
        # pnga4(paste(colnames(counts)[i], '-fit.png', sep=''))
        # for (j in 1:nrow(median.counts))
          # m[as.character(median.counts[j, 'mappability']), as.character(median.counts[j, 'gc'])] <- fit[j]
        # image(x, y, m, col=paste('#', c(sprintf('%02X', 0:255), rep('FF', 256)), c(rep('FF', 256), sprintf('%02X', 255:0)), sprintf('%02X', 255), sep=''), xlab='mappability', ylab='GC content', main=paste('Loess Fit -', colnames(counts)[i]), zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))))
        # contour(x, y, m, nlevels=20, zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))), add=TRUE)
        # dev.off()
        # pnga4(paste(colnames(counts)[i], '-residuals.png', sep=''))
        # for (j in 1:nrow(median.counts))
          # m[as.character(median.counts[j, 'mappability']), as.character(median.counts[j, 'gc'])] <- l$residuals[j]
        # image(x, y, m, col=paste('#', c(sprintf('%02X', 0:255), rep('FF', 256)), c(rep('FF', 256), sprintf('%02X', 255:0)), sprintf('%02X', 255), sep=''), xlab='mappability', ylab='GC content', main=paste('Loess Residuals -', colnames(counts)[i]), zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))))
        # contour(x, y, m, nlevels=20, zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))), add=TRUE)
        # dev.off()
      # }
    }, silent=TRUE)
    corrected[,i] <- corvals
  }
  corrected <- round(corrected, digits=2)
  # if (plotting) {
    # setwd('..')
    # pnga4(sub('-counts\\.(txt|RData)', '-residuals.png', basename(countsfile)))
    # image(x, y, m, col=paste('#', c(sprintf('%02X', 0:255), rep('FF', 256)), c(rep('FF', 256), sprintf('%02X', 255:0)), sprintf('%02X', 255), sep=''), xlab='mappability', ylab='GC content', main='Median Residuals', zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))))
    # contour(x, y, m, nlevels=20, zlim=c(-max(abs(range(m, finite=TRUE))), max(abs(range(m, finite=TRUE)))), add=TRUE)
    # dev.off()
  # }
  phenodata$loess.span <- used.span
  phenodata$loess.family <- used.family
  cat('Done.\n')
  list(phenodata=phenodata, bins=bins, counts=counts, corrected=corrected, residuals=residuals)
}

# a subfunction for performing the correction on one sample at a time could be defined here.

# EOF

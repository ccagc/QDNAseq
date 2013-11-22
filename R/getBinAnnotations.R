#########################################################################/**
# @RdocFunction getBinAnnotations
#
# @title "Gets bin annotation data for a particular bin size"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{binSize}{A @numeric scalar specifying the width of the bins
#    in units of kbp (1000 base pairs), e.g. \code{binSize=15} corresponds
#    to 15 kbp bins.}
#   \item{genome}{A @character string specify the genome and genome version
#    to be used.}
#   \item{type}{A @character string specify the experiment type, e.g. "SR50"
#    or "PE1000".}
#   \item{cache}{If @TRUE, the retrieved bin annotation data is cached
#    on the file system, otherwise not.}
#   \item{force}{If @TRUE, the bin anonnation data is retrieved/calculated
#    regardless of it already exists in the cache or not.}
# }
#
# \value{
#   Returns a @see "Biobase::AnnotatedDataFrame" object.
# }
#
# @author "IS"
#
# \seealso{
#   @see "createBins".
# }
#
# @keyword IO
#*/#########################################################################
getBinAnnotations <- function(binSize, genome='hg19', type='SR50', cache=TRUE,
  force=FALSE) {
  cacheKey <- list(genome=genome, binSize=binSize, type=type)
  cacheDir <- c('QDNAseq', 'binAnnotations')
  cacheSuffix <- paste('.', genome, '.', binSize, 'kbp.', type, sep='')
  if (!force) {
    bins <- loadCache(key=cacheKey, suffix=cacheSuffix, dirs=cacheDir)
    if (!is.null(bins)) {
      vmsg('Bin annotations for genome ', genome, ', bin size of ',
        binSize, 'kbp, and experiment type ', type, ' loaded from cache.')

      if (is.null(attr(bins, 'QDNAseqVersion')) ||
        attr(bins, 'QDNAseqVersion') < '0.5.4') {
        vmsg('Old version detected, ignoring.')
        bins <- NULL
      }
    }
  } else {
    bins <- NULL
  }


  if (is.null(bins)) {
    vmsg('Downloading bin annotations for genome ', genome,
      ', bin size ', binSize, 'kbp, and experiment type ', type, ' ...',
      appendLF=FALSE)
    filename <- sprintf('QDNAseq.%s.%gkbp.%s.rds', genome, binSize, type)
    urlPath <- 'http://cdn.bitbucket.org/ccagc/qdnaseq/downloads'
    remotefile <- file.path(urlPath, filename, fsep='/')
    localfile <- tempfile()
    tryCatch({
      result <- downloadFile(remotefile, localfile)
    }, error=function(e) {
      vmsg(' not found. Please generate them first.')
      stop(e)
    })
    bins <- readRDS(localfile)
    file.remove(localfile)
    if (cache) {
      vmsg(' saving in cache ...')
      saveCache(bins, key=cacheKey, suffix=cacheSuffix, dirs=cacheDir,
        compress=TRUE)
    }
    vmsg()
  }
  bins
}




#########################################################################/**
# @RdocFunction createBins
#
# @alias calculateMappability
# @alias calculateBlacklist
# @alias iterateResiduals
#
# @title "Builds bin annotation data for a particular bin size"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{bsgenome}{A BSgenome ...}
#   \item{binSize}{A @numeric scalar specifying the width of the bins
#    in units of kbp (1000 base pairs), e.g. \code{binSize=15} corresponds
#    to 15 kbp bins.}
#   \item{ignoreUnderscored}{Whether to ignore sequences with underscores
#     in their names ...}
#   \item{ignoreMitochondria}{Wheter to ignore mitochondrial DNA  ...}
# }
#
# \value{
#   Returns ...
# }
#
# @author "IS"
#
# \seealso{
#   @see "getBinAnnotations".
# }
#*/#########################################################################
createBins <- function(bsgenome, binSize, ignoreUnderscored=TRUE,
  ignoreMitochondria=TRUE) {
  chrs <- GenomicRanges::seqnames(bsgenome)
  if (ignoreUnderscored)
    chrs <- chrs[-grep('_', chrs)]
  if (ignoreMitochondria)
    chrs <- chrs[-grep('^(chr)?M(T)?$', chrs)]
  lengths <- GenomicRanges::seqlengths(bsgenome)[chrs]
  start <- end <- integer()
  bases <- gc <- numeric()
  vmsg('Creating bins of ', binSize, ' kbp for genome ',
    substitute(bsgenome))

  # Bin size in units of base pairs
  binWidth <- binSize*1000L

  for (chr in chrs) {
    vmsg('  Processing ', chr, ' ...', appendLF=FALSE)
    chr.size <- lengths[chr]
    chr.starts <- seq(from=1, to=chr.size, by=binWidth)
    chr.ends <- chr.starts + binWidth - 1L
    chr.ends[length(chr.ends)] <- chr.size
    chr.seq <- BSgenome::getSeq(bsgenome, chr, as.character=TRUE)
    bin.seq <- substring(chr.seq, first=chr.starts, last=chr.ends)
    acgt <- gsub('[^ACGT]', '', bin.seq)
    cg <- gsub('[^CG]', '', acgt)
    chr.bases <- nchar(acgt) / (binWidth) * 100
    chr.gc <- nchar(cg) / nchar(acgt) * 100
    start <- c(start, chr.starts)
    end <- c(end, chr.ends)
    bases <- c(bases, chr.bases)
    gc <- c(gc, chr.gc)
    vmsg()
  }
  gc[is.nan(gc)] <- NA_real_
  bins <- data.frame(chromosome=rep(chrs, times=ceiling(lengths/binWidth)),
                     start, end, bases, gc, stringsAsFactors=FALSE)
  bins$chromosome <- sub('^chr', '', bins$chromosome)
  rownames(bins) <- sprintf('%s:%i-%i', bins$chromosome, bins$start, bins$end)
  bins
}

calculateMappability <- function(bins, bigWigFile,
  bigWigAverageOverBed='bigWigAverageOverBed') {
  vmsg('Calculating mappabilities per bin from file\n  ', bigWigFile,
    '\n  ',
    appendLF=FALSE)
  binbed <- tempfile(fileext='.bed')
  mapbed <- tempfile(fileext='.bed')
  bins <- bins[,c('chromosome', 'start', 'end')]
  bins$chromosome <- paste('chr', bins$chromosome, sep='')
  bins$start <- bins$start - 1
  bins$name <- seq_len(nrow(bins))
  scipen <- options('scipen')
  options(scipen=10)
  write.table(bins, binbed, quote=FALSE, sep='\t', row.names=FALSE,
    col.names=FALSE)
  options(scipen=scipen)
  cmd <- paste(bigWigAverageOverBed, ' "', bigWigFile, '" "', binbed,
    '" -bedOut="', mapbed, '" /dev/null', sep='')
  system(cmd)
  map <- read.table(mapbed, sep='\t', as.is=TRUE)
  map$V5 * 100
}

calculateBlacklist <- function(bins, bedFiles, ncpus=1) {
  vmsg('Calculating overlaps per bin with BED files \n  ', paste(bedFiles,
    collapse='\n  '), '\n  ...', appendLF=FALSE)
  beds <- list()
  for (bed in bedFiles)
    beds[[bed]] <- read.table(bed, sep='\t', as.is=TRUE)
  combined <- beds[[1L]]
  if (length(beds) >= 2L)
    for (i in 2:length(beds))
      combined <- rbind(combined, beds[[i]])
  combined <- combined[, 1:3]
  colnames(combined) <- c('chromosome', 'start', 'end')
  combined$chromosome <- sub('^chr', '', combined$chromosome)
  combined <- combined[combined$chromosome %in% unique(bins$chromosome), ]
  combined$chromosome[combined$chromosome=='X'] <- '23'
  combined$chromosome[combined$chromosome=='Y'] <- '24'
  combined$chromosome <- as.integer(combined$chromosome)
  combined <- combined[!is.na(combined$chromosome), ]
  combined$start <- combined$start + 1
  combined <- combined[order(combined$chromosome, combined$start), ]
  joined <- data.frame()
  prev <- combined[1L,]
  # Sanity check
  stopifnot(nrow(combined) >= 2L);
  for (i in 2:nrow(combined)) {
    if (combined[i, 'chromosome'] != prev$chromosome ||
      combined[i, 'start'] > (prev$end + 1)) {
      joined <- rbind(joined, prev)
      prev <- combined[i,]
    } else {
      prev$end <- max(prev$end, combined[i, 'end'])
    }
  }
  joined <- rbind(joined, prev)
  bins$chromosome[bins$chromosome=='X'] <- '23'
  bins$chromosome[bins$chromosome=='Y'] <- '24'
  bins$chromosome <- as.integer(bins$chromosome)
  overlap.counter <- function(x, joined) {
    chr <- as.integer(x['chromosome'])
    start <- as.integer(x['start'])
    end <- as.integer(x['end'])
    overlaps <- joined[joined$chromosome   == chr &
                       joined$start        <= end &
                       joined$end          >= start, ]
    bases <- 0
    for (i in rownames(overlaps))
      bases <- bases + min(end, overlaps[i, 'end']) -
        max(start, overlaps[i, 'start']) + 1
    bases / (end - start + 1) * 100
  }
  if (ncpus > 1L) {
    snowfall::sfInit(parallel=TRUE, cpus=ncpus)
    snowfall::sfExport(list=c('overlap.counter'))
    blacklist <- snowfall::sfApply(bins, 1, overlap.counter, joined)
    snowfall::sfStop()
  } else {
    blacklist <- apply(bins, MARGIN=1L, FUN=overlap.counter, joined)
  }
  vmsg()
  blacklist
}

iterateResiduals <- function(object, cutoff=4.0, maxIter=30, ...) {
  first <- sum(binsToUse(object))
  previous <- first
  vmsg('Iteration #1 with ', format(previous, big.mark=','),
    ' bins.')
  object <- correctBins(object, storeResiduals=TRUE, ...)
  residuals <- assayDataElement(object, 'residuals')
  residuals[!binsToUse(object), ] <- NA
  residual <- apply(residuals, 1, median, na.rm=TRUE)
  cutoffValue <- cutoff * madDiff(residual, na.rm=TRUE)
  if (is.numeric(cutoff))
    binsToUse(object) <- binsToUse(object) & !is.na(residual) &
      abs(residual) <= cutoffValue
  num <- sum(binsToUse(object))
  iter <- 2
  while (previous != num && iter <= maxIter) {
    previous <- num
    vmsg('Iteration #', iter, ' with ', format(previous, big.mark=','),
      ' bins.')
    object <- correctBins(object, storeResiduals=TRUE, ...)
    residuals <- assayDataElement(object, 'residuals')
    residuals[!binsToUse(object), ] <- NA
    residual <- apply(residuals, 1, median, na.rm=TRUE)
    binsToUse(object) <- binsToUse(object) & !is.na(residual) &
      abs(residual) <= cutoffValue
    num <- sum(binsToUse(object))
    iter <- iter + 1
  }
  if (previous == num) {
    vmsg('Convergence at ', format(previous, big.mark=','), ' bins.')
    vmsg(format(first-previous, big.mark=','), ' additional bins removed.')
  } else if (iter == maxIter) {
    vmsg('Reached maxIter=', maxIter, ' iterations without convergence.')
  }
  residual
}

# EOF

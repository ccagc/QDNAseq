#########################################################################/**
# @RdocFunction binReadCounts
#
# @title "Calculate binned read counts from a set of BAM files"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{bins}{...}
#   \item{bamfiles}{...}
#   \item{path}{...}
#   \item{ext}{...}
#   \item{bamnames}{...}
#   \item{phenofile}{...}
#   \item{...}{Additional arguments passed to @see ".binReadCountsPerSample"}
# }
#
# \value{
#   Returns a named @list with elements ...
# }
#
# @author "IS"
#
# \seealso{
#   Internally, @see ".binReadCountsPerSample" is used.
# }
#
# @keyword IO
#*/#########################################################################
binReadCounts <- function(bins, bamfiles=NULL, path=NULL, ext='bam',
  bamnames=NULL, phenofile=NULL, ...) {
  if (is.null(bamfiles))
    bamfiles <- list.files(ifelse(is.null(path), '.', path),
      pattern=sprintf('%s$', ext), full.names=TRUE)
  if (length(bamfiles) == 0L)
    stop('No files to process.')
  if (is.null(bamnames)) {
    bamnames <- basename(bamfiles)
    bamnames <- sub(sprintf('[\\.]?%s$', ext), '', bamnames)
  } else if (length(bamfiles) != length(bamnames)) {
    stop('bamfiles and bamnames have to be of same length.')
  }
  phenodata <- data.frame(name=bamnames, row.names=bamnames,
    stringsAsFactors=FALSE)
  if (!is.null(phenofile)) {
    pdata <- read.table(phenofile, header=TRUE, sep='\t', as.is=TRUE,
      row.names=1L)
    phenodata <- cbind(phenodata, pdata[rownames(phenodata), ])
  }

  counts <- matrix(NA_integer_, nrow=nrow(bins), ncol=length(bamnames),
    dimnames=list(rownames(bins), bamnames))
  for (i in seq_along(bamfiles)) {
    counts[, i] <- .binReadCountsPerSample(bins=bins,
      bamfile=bamfiles[i], ...)
    gc(FALSE)
  }

  phenodata$reads <- colSums(counts)
  new('QDNAseqReadCounts', bins=bins, counts=counts, phenodata=phenodata)
}



#########################################################################/**
# @RdocFunction .binReadCountsPerSample
#
# @title "Calculate binned read counts from a BAM file"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{bins}{...}
#   \item{bamfile}{...}
#   \item{cache}{...}
#   \item{force}{...}
#   \item{maxChunk}{...}
#   \item{isPaired=NA}{...}
#   \item{isProperPair=NA}{...}
#   \item{isUnmappedQuery=FALSE}{...}
#   \item{hasUnmappedMate=NA}{...}
#   \item{isMinusStrand=NA}{...}
#   \item{isMateMinusStrand=NA}{...}
#   \item{isFirstMateRead=NA}{...}
#   \item{isSecondMateRead=NA}{...}
#   \item{isNotPrimaryRead=NA}{...}
#   \item{isNotPassingQualityControls=NA}{...}
#   \item{isDuplicate=FALSE}{...}
#   \item{minMapq}{If quality scores exists, the minimum quality score required in order to keep a read, otherwise all reads are kept.}
# }
#
# \value{
#   Returns an @integer @vector.
# }
#
# @author "IS"
#
# \seealso{
#   To retrieve bin read counts for a set of BAM files,
#   see @see "binReadCounts".
# }
#
# @keyword IO
# @keyword internal
#*/#########################################################################
.binReadCountsPerSample <- function(bins, bamfile, cache=TRUE, force=FALSE,
  maxChunk=100000000L, isPaired=NA, isProperPair=NA, isUnmappedQuery=FALSE,
  hasUnmappedMate=NA, isMinusStrand=NA, isMateMinusStrand=NA,
  isFirstMateRead=NA, isSecondMateRead=NA, isNotPrimaryRead=NA,
  isNotPassingQualityControls=FALSE, isDuplicate=FALSE, minMapq=37) {

  binSize <- (bins$end[1L]-bins$start[1L]+1)/1000

  bamfile <- normalizePath(bamfile)
  fullname <- sub('\\.[^.]*$', '', basename(bamfile))

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  readCountCacheKey <- list(bamfile=bamfile, isPaired=isPaired,
    isProperPair=isProperPair, isUnmappedQuery=isUnmappedQuery,
    hasUnmappedMate=hasUnmappedMate, isMinusStrand=isMinusStrand,
    isMateMinusStrand=isMateMinusStrand, isFirstMateRead=isFirstMateRead,
    isSecondMateRead=isSecondMateRead, isNotPrimaryRead=isNotPrimaryRead,
    isNotPassingQualityControls=isNotPassingQualityControls,
    isDuplicate=isDuplicate, minMapq=minMapq, binSize=binSize)
  readCountCacheDir <- c('QDNAseq', 'readCounts')
  readCountCacheSuffix <- paste('.', fullname, '.', binSize, 'kbp', sep='')
  if (!force) {
    readCounts <- loadCache(key=readCountCacheKey, sources=bamfile,
      suffix=readCountCacheSuffix, dirs=readCountCacheDir)
    if (!is.null(readCounts)) {
      message('Loaded binned read counts from cache for ', basename(bamfile),
        appendLF=FALSE)
      if (is.null(attr(readCounts, 'QDNAseqVersion'))) {
        attr(readCounts, 'QDNAseqVersion') <- packageVersion('QDNAseq')
        message(', re-caching with version number ...', appendLF=FALSE)
        saveCache(readCounts, key=readCountCacheKey, sources=bamfile,
          suffix=readCountCacheSuffix, dirs=readCountCacheDir, compress=TRUE)
      }
      message()
      return(readCounts)
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve counts per chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  readCacheKey <- list(bamfile=bamfile, isPaired=isPaired,
    isProperPair=isProperPair, isUnmappedQuery=isUnmappedQuery,
    hasUnmappedMate=hasUnmappedMate, isMinusStrand=isMinusStrand,
    isMateMinusStrand=isMateMinusStrand, isFirstMateRead=isFirstMateRead,
    isSecondMateRead=isSecondMateRead, isNotPrimaryRead=isNotPrimaryRead,
    isNotPassingQualityControls=isNotPassingQualityControls,
    isDuplicate=isDuplicate, minMapq=minMapq)
  readCacheDir <- c('QDNAseq', 'reads')
  readCacheSuffix <- paste('.', fullname, sep='')
  hits <- NULL
  if (!force)
    hits <- loadCache(key=readCacheKey, sources=bamfile,
      suffix=readCacheSuffix, dirs=readCacheDir)

  if (!is.null(hits)) {
    message('Loaded reads from cache for ', basename(bamfile), ',',
      appendLF=FALSE)
    if (is.null(attr(hits, 'QDNAseqVersion'))) {
      attr(hits, 'QDNAseqVersion') <- packageVersion('QDNAseq')
      message(' re-caching with version number ...', appendLF=FALSE)
      saveCache(hits, key=readCacheKey, sources=bamfile,
        suffix=readCacheSuffix, dirs=readCacheDir, compress=TRUE)
    }
  } else {
    message('Extracting reads from ', basename(bamfile), ' ...',
      appendLF=FALSE)
    flag <- scanBamFlag(isPaired=isPaired,
      isProperPair=isProperPair, isUnmappedQuery=isUnmappedQuery,
      hasUnmappedMate=hasUnmappedMate, isMinusStrand=isMinusStrand,
      isMateMinusStrand=isMateMinusStrand, isFirstMateRead=isFirstMateRead,
      isSecondMateRead=isSecondMateRead, isNotPrimaryRead=isNotPrimaryRead,
      isNotPassingQualityControls=isNotPassingQualityControls,
      isDuplicate=isDuplicate)
    params <- ScanBamParam(flag=flag, what=c('rname', 'pos', 'mapq'))
    reads <- scanBam(bamfile, param=params)
    reads <- reads[[1L]]

    # Filter by read quality scores?
    hasMapq <- any(is.finite(reads[['mapq']]))
    if (hasMapq) {
      keep <- which(reads[['mapq']] >= minMapq)
      reads <- lapply(reads, FUN=function(x) x[keep])
    }

    # Drop quality scores - not needed anymore
    reads[['mapq']] <- NULL

    # Sort counts by chromosome
    hits <- list()
    chrs <- unique(reads[['rname']])
    for (chr in chrs) {
      keep <- which(reads[['rname']] == chr)
      hits[[chr]] <- reads[['pos']][keep]
    }
    names(hits) <- sub('^chr', '', names(hits))
    rm(list=c('reads'))
    gc(FALSE)

    if (cache) {
      attr(hits, 'QDNAseqVersion') <- packageVersion('QDNAseq')
      message(' saving in cache ...', appendLF=FALSE)
      saveCache(hits, key=readCacheKey, sources=bamfile,
        suffix=readCacheSuffix, dirs=readCacheDir, compress=TRUE)
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Bin by chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  message(' binning ...', appendLF=FALSE)
  ## TO DO: the binning is very memory intensive, and therefore should be done
  ## to only a maximum of maxChunk reads at a time.
  readCounts <- integer(length=nrow(bins))
  for (chromosome in names(hits)) {
    keep <- which(bins$chromosome == chromosome);

    # No reads on this chromosome?
    if (length(keep) == 0L)
      next

    chromosome.breaks <- c(bins$start[keep], max(bins$end[keep]))
    ## without as.numeric(), command below gives an integer overflow warning:
    counts <- hist(hits[[chromosome]], breaks=as.numeric(chromosome.breaks),
      right=FALSE, plot=FALSE)$count

    ## count <- binCounts(hits[[chromosome]], chromosome.breaks) # matrixStats
    readCounts[keep] <- readCounts[keep] + counts

    ## Not needed anymore
    chromosome.breaks <- keep <- count <- NULL
  }
  ## Not needed anymore
  rm(list=c("hits"))
  gc(FALSE)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store results in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cache) {
    attr(readCounts, 'QDNAseqVersion') <- packageVersion('QDNAseq')
    message(' saving in cache ...', appendLF=FALSE)
    saveCache(readCounts, key=readCountCacheKey, sources=bamfile,
      suffix=readCountCacheSuffix, dirs=readCountCacheDir, compress=TRUE)
  }

  message()
  readCounts
}

# EOF

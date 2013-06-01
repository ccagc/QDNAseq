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
#   \item{allosomeBins}{...}
#   \item{incompleteBins}{...}
#   \item{blacklistedBins}{...}
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
binReadCounts <- function(bins, bamfiles=NULL, path='.', ext='bam', bamnames=NULL, phenofile=NULL, allosomeBins='flag', incompleteBins='flag', blacklistedBins='flag', ...) {
  if (is.null(bamfiles))
    bamfiles <- list.files(path, pattern=paste(ext, '$', sep=''))
  if (length(bamfiles) == 0L)
    stop('No files to process.')
  if (is.null(bamnames)) {
    bamnames <- sub(paste('\\.?', ext, '$', sep=''), '', bamfiles)
  } else if (length(bamfiles) != length(bamnames)) {
    stop('bamfiles and bamnames have to be of same length.')
  }
  phenodata <- data.frame(name=bamnames, row.names=bamnames, stringsAsFactors=FALSE)
  if (!is.null(phenofile)) {
    pdata <- read.table(phenofile, header=TRUE, sep='\t', as.is=TRUE, row.names=1L)
    phenodata <- cbind(phenodata, pdata[rownames(phenodata),])
  }
  counts <- matrix(nrow=nrow(bins), ncol=length(bamnames), dimnames=list(rownames(bins), bamnames))
  for (i in seq_along(bamfiles)) {
    counts[,i] <- .binReadCountsPerSample(bins=bins, bamfile=file.path(path, bamfiles[i]), ...)
    gc(FALSE)
  }
  phenodata$reads <- colSums(counts)
  condition <- condition <- rep(TRUE, times=nrow(bins))
  if (allosomeBins=='flag')
    condition <- condition & bins$chromosome %in% as.character(1:22)
  if (incompleteBins=='flag') {
    condition <- condition & bins$bases == 100
  } else if (incompleteBins=='adjust') {
    counts[bins$bases==0] <- NA # theoretically, BWA might place reads where reference is all Ns
    counts <- counts / bins$bases * 100
  }
  if (blacklistedBins=='flag')
    condition <- condition & bins$blacklist == 0
  bins$filter <- condition
  new('qdnaseq', bins=bins, counts=counts, phenodata=phenodata)
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
#   \item{minMapq=37}{...}
# }
#
# \value{
#   Returns ...
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
.binReadCountsPerSample <- function(bins, bamfile, cache=TRUE, maxChunk=100000000L, isPaired=NA, isProperPair=NA, isUnmappedQuery=FALSE, hasUnmappedMate=NA, isMinusStrand=NA, isMateMinusStrand=NA, isFirstMateRead=NA, isSecondMateRead=NA, isNotPrimaryRead=NA, isNotPassingQualityControls=FALSE, isDuplicate=FALSE, minMapq=37) {
  binsize <- (bins$end[1L]-bins$start[1L]+1)/1000
  linkTarget <- Sys.readlink(bamfile)
  if (linkTarget != '') {
    bamfile <- linkTarget
  }
  readCounts <- NULL
  if (cache==TRUE)
    readCounts <- loadCache(key=list(bamfile=bamfile, isPaired=isPaired, isProperPair=isProperPair, isUnmappedQuery=isUnmappedQuery, hasUnmappedMate=hasUnmappedMate, isMinusStrand=isMinusStrand, isMateMinusStrand=isMateMinusStrand, isFirstMateRead=isFirstMateRead, isSecondMateRead=isSecondMateRead, isNotPrimaryRead=isNotPrimaryRead, isNotPassingQualityControls=isNotPassingQualityControls, isDuplicate=isDuplicate, minMapq=minMapq, binsize=binsize), sources=bamfile, dirs='QDNAseq')
  if (!is.null(readCounts)) {
    cat('Loaded binned read counts from cache for ', basename(bamfile), '\n', sep='')
    return(readCounts)
  }
  hits <- NULL
  if (cache==TRUE)
    hits <- loadCache(key=list(bamfile=bamfile, isPaired=isPaired, isProperPair=isProperPair, isUnmappedQuery=isUnmappedQuery, hasUnmappedMate=hasUnmappedMate, isMinusStrand=isMinusStrand, isMateMinusStrand=isMateMinusStrand, isFirstMateRead=isFirstMateRead, isSecondMateRead=isSecondMateRead, isNotPrimaryRead=isNotPrimaryRead, isNotPassingQualityControls=isNotPassingQualityControls, isDuplicate=isDuplicate, minMapq=minMapq), sources=bamfile, dirs='QDNAseq')
  if (!is.null(hits)) {
    cat('Loaded reads from cache for ', basename(bamfile), ',', sep='')
  } else {
    cat('Extracting reads from ', basename(bamfile), ' ...', sep='')
    hitsfile <- tempfile()
    reads <- scanBam(bamfile, param=ScanBamParam(flag=scanBamFlag(isPaired=isPaired, isProperPair=isProperPair, isUnmappedQuery=isUnmappedQuery, hasUnmappedMate=hasUnmappedMate, isMinusStrand=isMinusStrand, isMateMinusStrand=isMateMinusStrand, isFirstMateRead=isFirstMateRead, isSecondMateRead=isSecondMateRead, isNotPrimaryRead=isNotPrimaryRead, isNotPassingQualityControls=isNotPassingQualityControls, isDuplicate=isDuplicate), what=c('rname', 'pos', 'mapq')))[[1]]
    hits <- list()
    for (chr in unique(reads[['rname']]))
      hits[[chr]] <- reads[['pos']][reads[['rname']]==chr & reads[['mapq']] >= minMapq]
    names(hits) <- sub('^chr', '', names(hits))
    rm(list=c('reads'))
    gc(FALSE)
    if (cache==TRUE || cache=='overwrite') {
      cat(' saving in cache ...', sep='')
      saveCache(hits, key=list(bamfile=bamfile, isPaired=isPaired, isProperPair=isProperPair, isUnmappedQuery=isUnmappedQuery, hasUnmappedMate=hasUnmappedMate, isMinusStrand=isMinusStrand, isMateMinusStrand=isMateMinusStrand, isFirstMateRead=isFirstMateRead, isSecondMateRead=isSecondMateRead, isNotPrimaryRead=isNotPrimaryRead, isNotPassingQualityControls=isNotPassingQualityControls, isDuplicate=isDuplicate, minMapq=minMapq), sources=bamfile, dirs='QDNAseq', compress=TRUE)
    }
  }
  cat(' binning ...', sep='')
  # TO DO: the binning is very memory intensive, and therefore should be done
  # to only a maximum of maxChunk reads at a time.
  readCounts <- numeric(length=nrow(bins))
  for (chromosome in names(hits)) {
    if (!chromosome %in% unique(bins$chromosome))
      next
    chromosome.breaks <- c(bins[bins$chromosome==chromosome, 'start'], max(bins[bins$chromosome==chromosome, 'end']))
    # without the as.numeric() the command below gives an integer overflow warning:
    count <- hist(hits[[chromosome]], breaks=as.numeric(chromosome.breaks), right=FALSE, plot=FALSE)$count
    # count <- binCounts(hits[[chromosome]], chromosome.breaks) # matrixStats
    readCounts[bins$chromosome==chromosome] <- readCounts[bins$chromosome==chromosome] + count
  }
  # Not needed anymore
  rm(list=c("hits"))
  gc(FALSE)
  if (cache==TRUE || cache=='overwrite') {
    cat(' saving in cache ...', sep='')
    saveCache(readCounts, key=list(bamfile=bamfile, isPaired=isPaired, isProperPair=isProperPair, isUnmappedQuery=isUnmappedQuery, hasUnmappedMate=hasUnmappedMate, isMinusStrand=isMinusStrand, isMateMinusStrand=isMateMinusStrand, isFirstMateRead=isFirstMateRead, isSecondMateRead=isSecondMateRead, isNotPrimaryRead=isNotPrimaryRead, isNotPassingQualityControls=isNotPassingQualityControls, isDuplicate=isDuplicate, minMapq=minMapq, binsize=binsize), sources=bamfile, dirs='QDNAseq', compress=TRUE)
  }
  cat('\n', sep='')
  readCounts
}

# EOF

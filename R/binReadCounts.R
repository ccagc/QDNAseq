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
    counts[,i] <- .binReadCountsPerSample(bins, bamfile=bamfiles[i], path=path, ...)
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
#   \item{path}{...}
#   \item{cache}{...}
#   \item{samtools}{...}
#   \item{f}{...}
#   \item{F}{...}
#   \item{q}{...}
#   \item{maxChunk}{...}
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
.binReadCountsPerSample <- function(bins, bamfile, path, cache=TRUE, samtools='samtools', f='', F='0x0404', q=37, maxChunk=100000000L) {
  binsize <- (bins$end[1L]-bins$start[1L]+1)/1000
  linkTarget <- Sys.readlink(file.path(path, bamfile))
  if (linkTarget != '') {
    bamfile <- basename(linkTarget)
    path <- dirname(linkTarget)
  }
  readCounts <- NULL
  if (cache)
    readCounts <- loadCache(key=list(bamfile=file.path(path, bamfile), f=f, F=F, q=q, binsize=binsize), dirs='QDNAseq')
  if (!is.null(readCounts)) {
cat('bins from cache\n')
    return(readCounts)
  }
  hits <- NULL
  hits <- loadCache(key=list(bamfile=file.path(path, bamfile), f=f, F=F, q=q), dirs='QDNAseq')
if (!is.null(hits))
cat('hits from cache\n')
  if (is.null(hits)) {
    hitsfile <- tempfile()
    # TO DO: This system call requires that 'samtools', 'cut', 'tr'
    # and 'gzip' are on the PATH and available.  It cannot be expected
    # that the latter three are available on Windows.
    # Alternatives (cross-platform):
    # (1) Rsamtools package on Bioconductor.
    # (2) The samtoolsView() function of future aroma.seq package
    system(paste(samtools, ' view -f "', f, '" -F "', F, '" -q ', q, ' "', file.path(path, bamfile), '" | cut -f3,4 | tr -d chr | gzip > "', hitsfile, '"', sep=''))
    reads <- as.data.frame(scan(hitsfile, what=list(chromosome=character(), pos=integer()), sep='\t', quiet=TRUE), stringsAsFactors=FALSE)
    hits <- list()
    for (chr in unique(reads$chromosome))
      hits[[chr]] <- reads[reads$chromosome==chr, 'pos']
    if (cache)
      saveCache(hits, key=list(bamfile=file.path(path, bamfile), f=f, F=F, q=q), dirs='QDNAseq', compress=TRUE)
  }
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
  rm(list=c("hits", "count"))
  gc(FALSE)
  if (cache)
    saveCache(readCounts, key=list(bamfile=file.path(path, bamfile), f=f, F=F, q=q, binsize=binsize), dirs='QDNAseq', compress=TRUE)
  readCounts
}

# EOF

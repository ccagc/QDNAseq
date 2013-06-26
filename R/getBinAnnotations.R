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
#   \item{binsize}{A @numeric scalar specifying ...}
#   \item{genome}{A @character string ...}
#   \item{cache}{A @logical ...}
#   \item{force}{A @logical ...}
# }
#
# \value{
#   Returns ...
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
getBinAnnotations <- function(binsize, genome='hg19', cache=TRUE,
  force=FALSE) {
  genome.build <- as.integer(gsub('[^0-9]', '', genome))
  if (genome.build %in% c(19, 37)) {
    genome.name <- 'hg19'
  } else if (genome.build %in% c(18, 36)) {
    genome.name <- 'hg18'
  } else {
    stop('Unknown genome: ', genome)
  }
  bins <- NULL
  cacheKey <- list(genome=genome.name, binsize=binsize)
  cacheDir <- c('QDNAseq', 'binAnnotations')
  cacheSuffix <- paste('.', genome.name, '.', binsize, 'kbp', sep='')
  if (!force)
    # TO DO: somehow check if file available online is newer than cached one?
    bins <- loadCache(key=cacheKey, suffix=cacheSuffix, dirs=cacheDir)
  if (!is.null(bins)) {
    message('Bin annotations for genome ', genome.name, ' and bin size of ',
      binsize, 'kbp loaded from cache.')
    return(bins)
  }
  message('Downloading bin annotations for genome ', genome.name,
    ' and bin size of ', binsize, 'kbp ...', appendLF=FALSE)
  remotefile <-
    paste('http://cdn.bitbucket.org/ccagc/qdnaseq/downloads/QDNAseq.',
    genome.name, '.', binsize, 'kbp.rds', sep='')
  localfile <- tempfile()
  tryCatch({
    library('R.utils')
    result <- downloadFile(remotefile, localfile)
  }, error=function(e) {
    message(' not found. Please generate them first.')
    stop(e)
  })
  bins <- readRDS(localfile)
  file.remove(localfile)
  if (cache) {
    message(' saving in cache ...')
    saveCache(bins, key=cacheKey, suffix=cacheSuffix, dirs=cacheDir,
      compress=TRUE)
  }
  message()
  bins
}




#########################################################################/**
# @RdocFunction createBins
#
# @alias calculateMappability
# @alias calculateBlacklist
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
#   \item{binsize}{A @numeric scalar specifying ...}
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
createBins <- function(bsgenome, binsize, ignoreUnderscored=TRUE,
  ignoreMitochondria=TRUE) {
  chrs <- GenomicRanges::seqnames(bsgenome)
  if (ignoreUnderscored)
    chrs <- chrs[-grep('_', chrs)]
  if (ignoreMitochondria)
    chrs <- chrs[-grep('^(chr)?M(T)?$', chrs)]
  lengths <- GenomicRanges::seqlengths(bsgenome)[chrs]
  start <- end <- integer()
  bases <- gc <- numeric()
  message('Creating bins of ', binsize, ' kbp for genome ',
    substitute(bsgenome))
  for (chr in chrs) {
    message('  Processing ', chr, ' ...', appendLF=FALSE)
    chr.size <- lengths[chr]
    chr.starts <- seq(from=1, to=chr.size, by=binsize*1000L)
    chr.ends <- chr.starts + binsize*1000L - 1L
    chr.ends[length(chr.ends)] <- chr.size
    chr.seq <- BSgenome::getSeq(bsgenome, chr, as.character=TRUE)
    bin.seq <- substring(chr.seq, first=chr.starts, last=chr.ends)
    acgt <- gsub('[^ACGT]', '', bin.seq)
    cg <- gsub('[^CG]', '', acgt)
    chr.bases <- nchar(acgt) / (binsize*1000L) * 100
    chr.gc <- nchar(cg) / nchar(acgt) * 100
    start <- c(start, chr.starts)
    end <- c(end, chr.ends)
    bases <- c(bases, chr.bases)
    gc <- c(gc, chr.gc)
    message()
  }
  gc[is.nan(gc)] <- NA
  bins <- data.frame(chromosome=rep(chrs, times=ceiling(lengths /
    (binsize*1000L))), start, end, bases, gc, stringsAsFactors=FALSE)
  bins$chromosome <- sub('^chr', '', bins$chromosome)
  rownames(bins) <- sprintf('%s:%i-%i', bins$chromosome, bins$start, bins$end)
  bins
}

calculateMappability <- function(bins, bigWigFile,
  bigWigAverageOverBed='bigWigAverageOverBed') {
  message('Calculating mappabilities per bin from file\n  ', bigWigFile,
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
  system(paste(bigWigAverageOverBed, ' "', bigWigFile, '" "', binbed,
    '" -bedOut="', mapbed, '" /dev/null', sep=''))
  map <- read.table(mapbed, sep='\t', as.is=TRUE)
  map$V5 * 100
}

calculateBlacklist <- function(bins, bedFiles, ncpus=1) {
  message('Calculating overlaps per bin with BED files \n  ', paste(bedFiles,
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
  message()
  blacklist
}

# 1000 Genomes residuals?

# EOF

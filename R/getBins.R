#########################################################################/**
# @RdocFunction getBins
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
getBins <- function(binsize, genome='hg19', cache=TRUE) {
  genome.build <- as.integer(gsub('[^0-9]', '', genome))
  if (genome.build %in% c(19, 37)) {
    genome.name <- 'hg19'
  } else if (genome.build %in% c(18, 36)) {
    genome.name <- 'hg18'
  } else {
    stop('Unknown genome: ', genome)
  }
  bins <- NULL
  if (cache==TRUE)
    # TO DO: somehow check if file available online is newer than cached one?
    bins <- loadCache(key=list(genome=genome, binsize=binsize), dirs='QDNAseq')
  if (!is.null(bins)) {
    cat('Bin annotations for genome ', genome.name, ' and bin size of ', binsize, 'kbp loaded from cache.\n', sep='')
    return(bins)
  }
  cat('Downloading bin annotations for genome ', genome.name, ' and bin size of ', binsize, 'kbp ...', sep='')
  remotefile <- paste('http://cdn.bitbucket.org/ccagc/qdnaseq/downloads/QDNAseq.', genome.name, '.', binsize, 'kbp.rds', sep='')
  localfile <- tempfile()
  error <- TRUE
  try({
    result <- downloadFile(remotefile, localfile)
    error <- FALSE
  }, silent=TRUE)
  if (error || is.null(result))
    stop('not found. Please generate them first.')
  bins <- readRDS(localfile)
  file.remove(localfile)
  if (cache==TRUE || cache=='overwrite') {
    cat(' saving in cache ...', sep='')
    saveCache(bins, key=list(genome=genome, binsize=binsize), dirs='QDNAseq', compress=TRUE)
  }
  cat('\n', sep='')
  bins
}




#########################################################################/**
# @RdocFunction createBins
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
#   \item{binsize}{A @numeric scalar specifying ...}
#   \item{genome}{A @character string ...}
# }
#
# \value{
#   Returns ...
# }
#
# @author "IS"
#
# \seealso{
#   @see "getBins".
# }
#*/#########################################################################
createBins <- function(binsize, genome='hg19') {
  genome.build <- as.integer(gsub('[^0-9]', '', genome))
  if (genome.build %in% c(19, 37)) {
    genome <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
    names(genome) <- c(1:22, 'X', 'Y')
  } else if (genome.build %in% c(18, 36)) {
    genome <- c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
    names(genome) <- c(1:22, 'X', 'Y')
  } else {
    stop('Unknown genome: ', genome)
  }
  start <- end <- integer()
  for (chromosome in names(genome)) {
    chromosome.size <- genome[chromosome]
    chromosome.starts <- seq(from=1, to=chromosome.size, by=binsize*1000)
    chromosome.ends <- chromosome.starts + binsize*1000 - 1
    chromosome.ends[length(chromosome.ends)] <- chromosome.size
    start <- c(start, chromosome.starts)
    end <- c(end, chromosome.ends)
  }
  bins <- data.frame(chromosome=rep(names(genome), times=ceiling(genome / (binsize*1000))), start, end, stringsAsFactors=FALSE)
  rownames(bins) <- paste(bins$chromosome, ':', bins$start, '-', bins$end, sep='')
  bins
}

# add functions for calculating for each bin:
# - GC content
# - mappability
# - overlap with ENCODE blacklist
# - 1000 Genomes residuals

# EOF

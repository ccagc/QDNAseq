binReadCounts <- function(binsize, files=NULL, path='.', ext='bam', names=NULL, genome='hg19', cache=TRUE, samtools='samtools', f='', F='0x0404', q=37, maxChunk=100000000) {
  genome.build <- as.integer(gsub('[^0-9]', '', genome))
  if (genome.build == 36 || genome.build == 18) {
    annfile <- file.path('/scratch/data/ref/hg18/MPScall', paste('QDNAseq.', binsize*1000, '.txt.gz', sep=''))
  } else {
    annfile <- file.path('/scratch/data/ref/hg19/MPScall', paste('QDNAseq.', binsize*1000, '.txt.gz', sep=''))
  }
  if (!file.exists(annfile))
    stop('No annotations file found for genome ', genome, ' and bin size of ', binsize, 'kbp. Please generate it first.')
  bins <- read.table(annfile, header=TRUE, sep='\t', as.is=TRUE)
  rownames(bins) <- paste(bins$chromosome, ':', bins$start, '-', bins$end, sep='')

  if (is.null(files))
    files <- list.files(path, paste('\\.', ext, '$', sep=''))
  if (is.null(names)) {
    names <- sub(paste('\\.', ext, '$', sep=''), '', files)
  } else if (length(files) != length(names)) {
    stop('files and names have to be of same length.')
  }

  counts <- matrix(nrow=nrow(bins), ncol=length(names), dimnames=list(rownames(bins), names))
  for (i in 1:length(files)) {
    counts[,i] <- .binReadCountsPerSample(binsize, files[i], path, bins, cache, samtools, f, F, q, maxChunk)
    gc(FALSE)
  }
  dat <- list(bins=bins, counts=counts) # phenodata?
}

.binReadCountsPerSample <- function(binsize, bamfile, path, bins, cache, samtools, f, F, q, maxChunk) {
  linkTarget <- Sys.readlink(file.path(path, bamfile))
  if (linkTarget != '') {
    bamfile <- basename(linkTarget)
    path <- dirname(linkTarget)
  }
  if (cache & !file.exists(file.path(path, '.QDNAseq')))
    dir.create(file.path(path, '.QDNAseq'))
  binfile <- file.path(path, '.QDNAseq', paste(bamfile, '.f', f, '.F', F, '.q', q, '.', binsize*1000, '.txt.gz', sep=''))
  if (file.exists(binfile)) {
    readCounts <- scan(binfile, what=integer(), quiet=TRUE)
  } else {
    if (cache) {
      hitsfile <- file.path(path, '.QDNAseq', paste(bamfile, '.f', f, '.F', F, '.q', q, '.txt.gz', sep=''))
    } else {
      hitsfile <- tempfile()
    }
    if (!file.exists(hitsfile))
      system(paste(samtools, ' view -f "', f, '" -F "', F, '" -q ', q, ' "', file.path(path, bamfile), '" | cut -f3,4 | tr -d chr | gzip > "', hitsfile, '"', sep=''))
    chromosome.sizes <- aggregate(bins$end, by=list(chromosome=bins$chromosome), max)
    rownames(chromosome.sizes) <- chromosome.sizes$chromosome
    readCounts <- numeric(length=nrow(bins))
    skip <- 0
    while(1) {
      hits <- as.data.frame(scan(hitsfile, what=list(chromosome=character(), pos=integer()), sep='\t', nmax=maxChunk, skip=skip, quiet=TRUE), stringsAsFactors=FALSE)
      if (nrow(hits) == 0)
        break
      for (chromosome in unique(hits$chromosome)) {
        if (!chromosome %in% unique(bins$chromosome))
          next
        chromosome.breaks <- c(bins[bins$chromosome==chromosome, 'start'], chromosome.sizes[chromosome, 'x'])
        count <- hist(hits[hits$chromosome==chromosome, 'pos'], breaks=as.numeric(chromosome.breaks), plot=FALSE)$count # without the as.numeric() gives an integer overflow warning
        readCounts[bins$chromosome==chromosome] <- readCounts[bins$chromosome==chromosome] + count
      }
      skip <- skip + maxChunk
      rm(hits, count); gc(FALSE)
    }
    if (cache) {
      cat(readCounts, file=sub('\\.gz$', '', binfile), sep='\n')
      system(paste('gzip "', sub('\\.gz$', '', binfile), '"', sep=''))
    }
  }
  readCounts
}

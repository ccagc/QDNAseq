options(scipen=10) # change filenames to 'kbp'

binReadCounts <- function(bins, bamfiles=NULL, path='.', ext='bam', bamnames=NULL, phenofile=NULL, genome='hg19', cache=TRUE, samtools='samtools', f='', F='0x0404', q=37, maxChunk=100000000) {
  if (is.null(bamfiles))
    bamfiles <- list.files(path, paste('\\.', ext, '$', sep=''))
  if (length(bamfiles) == 0)
    stop('No files to process.')
  if (is.null(bamnames)) {
    bamnames <- sub(paste('\\.', ext, '$', sep=''), '', bamfiles)
  } else if (length(bamfiles) != length(bamnames)) {
    stop('bamfiles and bamnames have to be of same length.')
  }
  phenodata <- data.frame(name=bamnames, row.names=bamnames, stringsAsFactors=FALSE)
  if (!is.null(phenofile)) {
    pdata <- read.table(phenofile, header=TRUE, sep='\t', as.is=TRUE, row.names=1)
    phenodata <- cbind(phenodata, pdata[rownames(phenodata),])
  }
  counts <- matrix(nrow=nrow(bins), ncol=length(bamnames), dimnames=list(rownames(bins), bamnames))
  for (i in 1:length(bamfiles)) {
    counts[,i] <- .binReadCountsPerSample(bins, bamfiles[i], path, cache, samtools, f, F, q, maxChunk)
    gc(FALSE)
  }
  phenodata$reads <- apply(counts, 2, sum)
  dat <- list(phenodata=phenodata, bins=bins, counts=counts)
}

.binReadCountsPerSample <- function(bins, bamfile, path, cache, samtools, f, F, q, maxChunk) {
  binsize <- (bins$end[1]-bins$start[1]+1)/1000
  linkTarget <- Sys.readlink(file.path(path, bamfile))
  if (linkTarget != '') {
    bamfile <- basename(linkTarget)
    path <- dirname(linkTarget)
  }
  if (cache & !file.exists(file.path(path, '.QDNAseq')))
    dir.create(file.path(path, '.QDNAseq'))
  binfile <- file.path(path, '.QDNAseq', paste(bamfile, '.f', f, '.F', F, '.q', q, '.', binsize, 'kbp.txt.gz', sep=''))
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
    readCounts <- numeric(length=nrow(bins))
    skip <- 0
    while(1) {
      hits <- as.data.frame(scan(hitsfile, what=list(chromosome=character(), pos=integer()), sep='\t', nmax=maxChunk, skip=skip, quiet=TRUE), stringsAsFactors=FALSE)
      if (nrow(hits) == 0)
        break
      for (chromosome in unique(hits$chromosome)) {
        if (!chromosome %in% unique(bins$chromosome))
          next
        chromosome.breaks <- c(bins[bins$chromosome==chromosome, 'start'], max(bins[bins$chromosome==chromosome, 'end']))
        # without the as.numeric() the command below gives an integer overflow warning:
        count <- hist(hits[hits$chromosome==chromosome, 'pos'], breaks=as.numeric(chromosome.breaks), plot=FALSE)$count
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

# EOF

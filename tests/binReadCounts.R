library(QDNAseq)

if (requireNamespace("QDNAseq.hg19", quietly = TRUE)) {
  bins <- getBinAnnotations(500, genome = "hg19")
  print(bins)

  bam <- system.file("extdata", "ex1.bam", package = "Rsamtools")
  print(bam)

  counts <- binReadCounts(bins, bamfiles = bam)
  print(counts)

  ## BUG: https://github.com/ccagc/QDNAseq/issues/89
  ## counts2 <- binReadCounts(bins, bamfiles = bam, chunkSize = 10e3)
  ## print(counts2)
}



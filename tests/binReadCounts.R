library(QDNAseq)

if (requireNamespace("QDNAseq", quietly = TRUE)) {
  bins <- getBinAnnotations(500)
  print(bins)

  bamFile <- system.file("extdata", "ex1.bam", package="Rsamtools")
  print(bamFile)

  counts <- binReadCounts(bins, bamfiles = bamFile)
  print(counts)

  ## BUG: https://github.com/ccagc/QDNAseq/issues/89
  ## counts2 <- binReadCounts(bins, bamfiles = bamFile, chunkSize = 10e3)
  ## print(counts2)
}



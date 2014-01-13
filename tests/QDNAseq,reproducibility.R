######################################################################
# This scripts asserts that for each processing step of QDNAseq
# the output/results are reproducible (numerically equal).
######################################################################
library("QDNAseq")

# Load data
data(LGG150)
data <- LGG150

# Filter out "bad" bins
dataF <- applyFilters(data, residual=TRUE, blacklist=TRUE)
dataFr <- applyFilters(data, residual=TRUE, blacklist=TRUE)
stopifnot(all.equal(dataFr, dataF))

# Correct read counts as a function of GC content and mappability
dataC <- correctBins(dataF)
dataCr <- correctBins(dataF)
stopifnot(all.equal(dataCr, dataC))

# Normalize binned read counts to have diploid normal copy number
dataN <- normalizeBins(dataC)
dataNr <- normalizeBins(dataC)
stopifnot(all.equal(dataNr, dataN))

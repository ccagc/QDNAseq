######################################################################
# This scripts asserts that for each processing step of QDNAseq
# the output/results are reproducible (numerically equal).
######################################################################
library("QDNAseq")
options("QDNAseq::verbose"=FALSE)

# Load data
data(LGG150)
data <- LGG150

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TRUTH
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Filter out "bad" bins
dataF <- applyFilters(data, residual=TRUE, blacklist=TRUE)

# Correct read counts as a function of GC content and mappability
dataC <- correctBins(dataF)

# Normalize binned read counts to have diploid normal copy number
dataN <- normalizeBins(dataC)

# Segment copy numbers
set.seed(42)  ## segmentBins() relies on RNG via DNAcopy::segment()
fit <- segmentBins(dataN)

# Call copy-number segments
fitC <- callBins(fit)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# REPRODUCIBILITY
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
strategies <- c("sequential", "multisession")
if (parallelly::supportsMulticore()) strategies <- c(strategies, "multicore")

oplan <- future::plan("list")
for (strategy in strategies) {
  message(sprintf("Reproducibility with plan(\"%s\") ...", strategy))
  
  future::plan(strategy)
  
  dataFr <- applyFilters(data, residual=TRUE, blacklist=TRUE)
  stopifnot(all.equal(dataFr, dataF))
  
  dataCr <- correctBins(dataF)
  stopifnot(all.equal(dataCr, dataC))
  
  dataNr <- normalizeBins(dataC)
  stopifnot(all.equal(dataNr, dataN))
  
  set.seed(42)  ## segmentBins() relies on RNG via DNAcopy::segment()
  fitr <- segmentBins(dataNr)
  stopifnot(all.equal(fitr, fit))
  
  fitCr <- callBins(fitr)
  stopifnot(all.equal(fitCr, fitC))
  
  message(sprintf("Reproducibility with plan(\"%s\") ... done", strategy))
}

future::plan(oplan)

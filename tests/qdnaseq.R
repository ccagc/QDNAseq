library("qdnaseq")

# Load data
data(LGG150)
data <- LGG150
print(data)

# Correct read counts as a function of GC content and mappability
readCountPlot(data)
dataC <- correctBins(data)
print(dataC)

# Normalize binned read counts to have diploid normal copy number
dataN <- normalizeBins(dataC)
print(dataN)
plot(dataN)
highlightFilters(dataN, mappability=50, blacklist=0, residual=2, bases=0)

# Filter out "bad" bins
dataNF <- applyFilters(dataN, mappability=50, blacklist=0, residual=2, bases=0)
print(dataNF)

# Segment copy numbers
fit <- segmentBins(dataNF)
print(fit)
plot(fit)

# Call copy-number segments
fitC <- callBins(fit)
print(fitC)
plot(fitC)

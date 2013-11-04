library("QDNAseq")

# Load data
data(LGG150)
data <- LGG150
print(data)

# Plot isobars of read counts
readCountPlot(data)

# Plot copy number profile
plot(data, ylim=c(-100, 200))
highlightFilters(data, residual=TRUE, blacklist=TRUE)

# Filter out "bad" bins
dataF <- applyFilters(data, residual=TRUE, blacklist=TRUE)
print(dataF)
plot(dataF, ylim=c(-100, 200))

# Correct read counts as a function of GC content and mappability
dataC <- correctBins(dataF)
print(dataC)
plot(dataC, ylim=c(-100, 200))

# Normalize binned read counts to have diploid normal copy number
dataN <- normalizeBins(dataC)
print(dataN)
plot(dataN)

# Plot noise
noisePlot(dataN)

# Segment copy numbers
fit <- segmentBins(dataN)
print(fit)
plot(fit)

# Call copy-number segments
#fitC <- callBins(fit)
#print(fitC)
#plot(fitC)

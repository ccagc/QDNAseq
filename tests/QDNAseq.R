library("QDNAseq")

# Load data
data(LGG150)
data <- LGG150
print(data)
stopifnot(inherits(data, "QDNAseqReadCounts"))

# Plot isobars of read counts
readCountPlot(data)

# Plot copy number profile
plot(data, ylim=c(-100, 200))
highlightFilters(data, residual=TRUE, blacklist=TRUE)

# Filter out "bad" bins
dataF <- applyFilters(data, residual=TRUE, blacklist=TRUE)
print(dataF)
plot(dataF, ylim=c(-100, 200))
stopifnot(inherits(dataF, "QDNAseqReadCounts"))

# Correct read counts as a function of GC content and mappability
dataC <- correctBins(dataF)
print(dataC)
plot(dataC, ylim=c(-100, 200))
stopifnot(inherits(dataC, "QDNAseqCopyNumbers"))

# Normalize binned read counts to have diploid normal copy number
dataN <- normalizeBins(dataC)
print(dataN)
plot(dataN)
stopifnot(inherits(dataN, "QDNAseqCopyNumbers"))

# Plot noise
noisePlot(dataF)

# Segment copy numbers
fit <- segmentBins(dataN)
print(fit)
plot(fit)
stopifnot(inherits(fit, "QDNAseqCopyNumbers"))

# Call copy-number segments
#fitC <- callBins(fit)
#print(fitC)
#plot(fitC)

library(QDNAseq)
library(utils)

# Load data
data(LGG150)
data <- LGG150
print(data)
stopifnot(inherits(data, "QDNAseqReadCounts"))

# Exporting
# FIXME:
# formats = c("vcf", "seg") require 'calls'; give informative error message
formats <- c("tsv", "igv", "bed")
for (format in formats) {
  file <- tempfile(fileext = sprintf(".%s", format))
  exportBins(data, file = file, format = format)
  stopifnot(file_test("-f", file))
  file.remove(file)
}

# Plot isobars of read counts
isobarPlot(data)

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

# Exporting
# FIXME:
# formats = c("vcf", "seg") require 'calls'; give informative error message
formats <- c("tsv", "igv", "bed")
for (format in formats) {
  file <- tempfile(fileext = sprintf(".%s", format))
  exportBins(dataC, file = file, format = format)
  stopifnot(file_test("-f", file))
  file.remove(file)
}

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

# Exporting
# FIXME:
# formats = c("vcf", "seg") require 'calls'; give informative error message
formats <- c("tsv", "igv", "bed")
for (format in formats) {
  file <- tempfile(fileext = sprintf(".%s", format))
  exportBins(fit, file = file, format = format)
  stopifnot(file_test("-f", file))
  file.remove(file)
}

# Call copy-number segments
fitC <- callBins(fit)
print(fitC)
plot(fitC)

# Exporting
formats <- c("tsv", "igv", "bed", "vcf", "seg")
for (format in formats) {
  file <- tempfile(fileext = sprintf(".%s", format))
  exportBins(fitC, file = file, format = format)
  ## FIXME: formats = c("vcf", "seg") does not respect 'file'
  if (!format %in% c("seg", "vcf")) {
    stopifnot(file_test("-f", file))
    file.remove(file)
  }
}

library(QDNAseq)
library(Biobase) ## combine(), sampleNames()
library(utils)

do_plot <- TRUE
do_cleanup <- TRUE

# Load data
data(LGG150)
data <- LGG150
print(data)
stopifnot(inherits(data, "QDNAseqReadCounts"))

# Plot isobars of read counts
if (do_plot) isobarPlot(data)

# Plot copy number profile
if (do_plot) {
  plot(data, ylim=c(-100, 200))
  highlightFilters(data, residual=TRUE, blacklist=TRUE)
}

# Filter out "bad" bins
dataF <- applyFilters(data, residual=TRUE, blacklist=TRUE)
print(dataF)
if (do_plot) plot(dataF, ylim=c(-100, 200))
stopifnot(inherits(dataF, "QDNAseqReadCounts"))

# Correct read counts as a function of GC content and mappability
dataC <- correctBins(dataF)
print(dataC)
if (do_plot) plot(dataC, ylim=c(-100, 200))
stopifnot(inherits(dataC, "QDNAseqCopyNumbers"))

# Normalize binned read counts to have diploid normal copy number
dataN <- normalizeBins(dataC)
print(dataN)
if (do_plot) plot(dataN)
stopifnot(inherits(dataN, "QDNAseqCopyNumbers"))

# Plot noise
if (do_plot) noisePlot(dataF)

# Segment copy numbers
fit <- segmentBins(dataN)
print(fit)
if (do_plot) plot(fit)
stopifnot(inherits(fit, "QDNAseqCopyNumbers"))

# Call copy-number segments
fitC <- callBins(fit)
print(fitC)
if (do_plot) plot(fitC)


# ---------------------------------------------------------------
# Exporting
# ---------------------------------------------------------------
message("* exportBins() ...")

sets <- list(data = data, dataC = dataC, fit = fit, fitC = fitC)
for (name in names(sets)) {
  set <- sets[[name]]
  formats <- c("tsv", "igv", "bed")
  if (name == "fitC") formats <- c(formats, "vcf", "seg")
  for (format in formats) {
    types <- c("copynumber")
    if (name %in% c("fit", "fitC")) types <- c(types, "segments")
    if (name == "fitC") types <- c(types, "calls")
    for (type in types) {
      fileext <- sprintf(".%s.%s", type, format)
      templates <- c("QDNAseq-%s.", "QDNAseq-%i.", "QDNAseq-%03i.")
      if (ncol(set) == 1L || !(format %in% c("bed", "seg", "vcf"))) {
        templates <- c("QDNAseq.", templates)
      }
      for (template in templates) {
        file <- tempfile(pattern = template, fileext = fileext)
        message(sprintf("  - exportBins(<%d samples>, format=\"%s\", type=\"%s\", file=\"%s\")", ncol(set), format, type, template))
        file <- exportBins(set, file = file, format = format, type = type)
        message(sprintf("    File(s) written: [n=%d] %s",
                length(file), paste(sQuote(file), collapse = ", ")))
        stopifnot(all(file_test("-f", file)))
        if (do_cleanup) {
	  file.remove(file)
          stopifnot(!any(file_test("-f", file)))
	}
      }
    }
  }
}

sets <- list(data = data, dataC = dataC, fit = fit, fitC = fitC)
sets <- lapply(sets, FUN = function(set) {
  stopifnot(ncol(set) == 1L)
  name <- sampleNames(set)
  setA <- set
  sampleNames(setA) <- sprintf("%sa", name)
  setB <- set
  sampleNames(setB) <- sprintf("%sb", name)
  combine(setA, setB)
})

for (name in names(sets)) {
  set <- sets[[name]]
  stopifnot(ncol(set) == 2L)
  formats <- c("tsv", "igv", "bed")
  if (name == "fitC") formats <- c(formats, "vcf", "seg")
  for (format in formats) {
    types <- c("copynumber")
    if (name %in% c("fit", "fitC")) types <- c(types, "segments")
    if (name == "fitC") types <- c(types, "calls")
    for (type in types) {
      fileext <- sprintf(".%s.%s", type, format)
      templates <- c("QDNAseq-%s.", "QDNAseq-%i.", "QDNAseq-%03i.")
      if (ncol(set) == 1L || !(format %in% c("bed", "seg", "vcf"))) {
        templates <- c("QDNAseq.", templates)
      }
      for (template in templates) {
        file <- tempfile(pattern = template, fileext = fileext)
        message(sprintf("  - exportBins(<%d samples>, format=\"%s\", type=\"%s\", file=\"%s\")", ncol(set), format, type, template))
        file <- exportBins(set, file = file, format = format, type = type)
        message(sprintf("    File(s) written: [n=%d] %s",
                length(file), paste(sQuote(file), collapse = ", ")))
        stopifnot(all(file_test("-f", file)))
        if (do_cleanup) {
          file.remove(file)
          stopifnot(!any(file_test("-f", file)))
	}  
      }
    }
  }
}

message("* exportBins() ... done")

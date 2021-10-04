library(QDNAseq)
library(Biobase)
library(utils)
set.seed(0xBEEF)

data(LGG150)
data <- LGG150
print(data)

dataF <- applyFilters(data, residual=TRUE, blacklist=TRUE)
print(dataF)

dataC <- correctBins(dataF)

## Force results to have segments gain, copy neutral, and gain.
cnAll <- assayDataElement(dataC, "copynumber")
cnAll[,1] <- rnorm(nrow(cnAll), mean = 1.0, sd = 0.05)
chr7 <- (chromosomes(dataC) == "7")
cn <- cnAll[chr7, , drop = FALSE]
n <- nrow(cn)
idxs <- seq(from=1/3*n - 0.1*n, to=1/3*n + 0.1*n)
cn[idxs,1] <- rnorm(length(idxs), mean = 2.0, sd = 0.05)
idxs <- seq(from=2/3*n - 0.1*n, to=2/3*n + 0.1*n)
cn[idxs,1] <- rnorm(length(idxs), mean = 2.0, sd = 0.05)
cnAll[chr7, ] <- cn
assayDataElement(dataC, "copynumber") <- cnAll
print(dataC)

dataN <- normalizeBins(dataC)
print(dataN)

fit <- segmentBins(dataN)
print(fit)

fitC <- callBins(fit)
print(fitC)

## Assert that everything is called copy-neutral
calls <- assayDataElement(fitC, "calls")
stopifnot(all(is.na(calls) | calls %in% c(0, 2)))

formats <- c("tsv", "igv", "bed", "seg", "vcf")
types <- c("copynumber", "segments", "calls")
for (format in formats) {
  for (type in types) {
    fileext <- sprintf(".%s.%s", type, format)
    file <- tempfile(pattern =  "QDNAseq-%s", fileext = fileext)
    message(sprintf("  - exportBins(<%d samples>, format=\"%s\", type=\"%s\")", ncol(fitC), format, type))
    file <- exportBins(fitC, format = format, type = type, file = file)
    message(sprintf("    File(s) written: [n=%d] %s",
            length(file), paste(sQuote(file), collapse = ", ")))
    stopifnot(all(file_test("-f", file)))

    if (format == "seg") {
      segs <- read.table(file, sep = "\t", header = TRUE)
      print(segs)
      stopifnot(all(segs$CHROMOSOME == "7"), nrow(segs) == 2L)
    } else if (format == "vcf") {
      segs <- read.table(file, sep = "\t", header = FALSE)
      print(segs)
      stopifnot(all(segs$V1 == "7"), nrow(segs) == 2L)
    }

    file.remove(file)
    stopifnot(!any(file_test("-f", file)))
  }
}

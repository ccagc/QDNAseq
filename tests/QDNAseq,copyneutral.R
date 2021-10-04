library(QDNAseq)
library(Biobase)
set.seed(0xBEEF)

data(LGG150)
data <- LGG150
print(data)

dataF <- applyFilters(data, residual=TRUE, blacklist=TRUE)
print(dataF)

dataC <- correctBins(dataF)

## Force all copy neutral data
cn <- assayDataElement(dataC, "copynumber")
cn[,1] <- rnorm(nrow(cn), mean = 1.0, sd = 0.05)
assayDataElement(dataC, "copynumber") <- cn

print(dataC)

dataN <- normalizeBins(dataC)
print(dataN)

fit <- segmentBins(dataN)
print(fit)

fitC <- callBins(fit)
print(fitC)


## Assert that everything is called copy-neutral
calls <- assayDataElement(fitC, "calls")
stopifnot(all(is.na(calls) | calls == 0))

formats <- c("tsv", "igv", "bed", "seg", "vcf")
types <- c("copynumber", "segments", "calls")
for (format in formats) {
  for (type in types) {
    fileext <- sprintf(".%s.%s", type, format)
    file <- tempfile(pattern =  "QDNAseq-%s", fileext = fileext)
    message(sprintf("  - exportBins(<%d samples>, format=\"%s\", type=\"%s\")", ncol(fitC), format, type))
    exportBins(fitC, format = format, type = type, file = file)
    message(sprintf("    File(s) written: [n=%d] %s",
            length(file), paste(sQuote(file), collapse = ", ")))
    stopifnot(all(file_test("-f", file)))
    file.remove(file)
    stopifnot(!any(file_test("-f", file)))
  }
}

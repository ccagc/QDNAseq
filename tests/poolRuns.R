library(QDNAseq)
library(Biobase)  ## sampleNames()

data(LGG150)

message("poolRuns(LGG150) ...")
pooledReadCounts <- poolRuns(LGG150, samples = sampleNames(LGG150))
print(pooledReadCounts)
message("poolRuns(LGG150) ... done")


message("Create fake data set with two samples ...")
phenodata <- phenoData(LGG150)
phenodata@data <- rbind(phenodata@data, phenodata@data)

counts <- assayDataElement(LGG150, "counts")
counts2 <- cbind(counts, counts)
colnames(counts2) <- sampleNames(phenodata)

x <- new("QDNAseqReadCounts", bins = featureData(LGG150), phenodata = phenodata, counts = counts2)
print(x)
message("Create fake data set with two samples ... done")

message("poolRuns(LGG150set) ...")
x_pool <- poolRuns(x, samples = sampleNames(x))
print(x_pool)

if (FALSE) {
message("poolRuns(LGG150set) - pooled ...")
## BUG #112 (https://github.com/ccagc/QDNAseq/issues/112)
x_pool <- poolRuns(x, samples = rep("pooled", times = 2L))
print(x_pool)
}

message("poolRuns(LGG150set) ... done")

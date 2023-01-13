library(QDNAseq)

data(LGG150)
readCounts <- LGG150

# Note: the following command will "pool" data from a single run, which
# does not really test much
pooledReadCounts <- poolRuns(readCounts, samples = "LGG150")
print(pooledReadCounts)

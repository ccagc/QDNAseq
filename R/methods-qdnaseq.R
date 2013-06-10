setMethod("initialize", "qdnaseq", function(.Object, bins, counts, phenodata,
  ...) {
  callNextMethod(.Object, featureData=AnnotatedDataFrame(bins),
    assayData=assayDataNew(counts=counts),
    phenoData=AnnotatedDataFrame(phenodata), ...)
})

# EOF

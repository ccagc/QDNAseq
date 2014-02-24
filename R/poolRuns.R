#########################################################################/**
# @RdocFunction poolRuns
#
# @alias poolRuns,QDNAseqSignals,character-method
#
# @title "Pools binned read counts across samples"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{object}{...}
#   \item{samples}{...}
#   \item{force}{...}
# }
#
# \value{
#   Returns a named @list containing elements ...
# }
#
# @author "IS"
#
# \seealso{
#   Internally, ...
# }
#
#*/#########################################################################
setMethod("poolRuns", signature=c(object="QDNAseqSignals",
  samples="character"), definition=function(object, samples, force=FALSE) {

  if (!force && "segmented" %in% assayDataElementNames(object))
    stop("Data has already been segmented. Pooling runs will ",
      "remove segmentation (and possible calling) ",
      "results. Please specify force=TRUE, if you want this.")

  if (length(samples) != ncol(object))
    stop("Parameter samples must be a vector of equal length as there are ",
      " samples in object.")

  newsamples <- sort(unique(samples))
  if (length(newsamples) == ncol(object)) {
    sampleNames(object) <- samples
    object <- object[,order(samples)]
    return(object)
  }
  oldsamples <- sampleNames(object)

  bins <- featureData(object)
  phenodata <- pData(object)
  phenodata[, 1] <- samples
  newphenodata <- phenodata[0, ]
  if (inherits(object, "QDNAseqReadCounts")) {
    counts <- assayDataElement(object, "counts")
    newcounts <- matrix(NA_integer_, nrow=nrow(object),
      ncol=length(newsamples), dimnames=list(featureNames(object), newsamples))
  } else if (inherits(object, "QDNAseqCopyNumbers")) {
    copynumber <- assayDataElement(object, "copynumber")
    newcopynumber <- matrix(NA_real_, nrow=nrow(object),
      ncol=length(newsamples), dimnames=list(featureNames(object), newsamples))
  }
  
  concatenateIfNotEqual <- function(x) {
    x <- sort(unique(x))
    paste(x, collapse=";")
  }

  for (newsample in newsamples) {
    replicates <- samples == newsample
    if (inherits(object, "QDNAseqReadCounts")) {
      newcounts[, newsample] <- rowSums(counts[, replicates, drop=FALSE])
    } else if (inherits(object, "QDNAseqCopyNumbers")) {
      newcopynumber[, newsample] <- rowMeans(copynumber[, replicates,
        drop=FALSE])
    }
    oldphenodata <- phenodata[replicates, ]
    totalReads <- sum(oldphenodata$total.reads)
    usedReads <- sum(oldphenodata$used.reads)
    numericCols <- sapply(oldphenodata, is.numeric)
    oldphenodata[1, numericCols] <- apply(oldphenodata[, numericCols,
      drop=FALSE], 2, mean)
    oldphenodata[1, !numericCols] <- apply(oldphenodata[, !numericCols,
      drop=FALSE], 2, concatenateIfNotEqual)
    oldphenodata[1, "total.reads"] <- totalReads
    oldphenodata[1, "used.reads"] <- usedReads
    oldphenodata[1, "expected.variance"] <- sum(binsToUse(object)) / usedReads
    
    newphenodata <- rbind(newphenodata, oldphenodata[1,])
  }
  rownames(newphenodata) <- newphenodata[, 1]
  newphenodata <- AnnotatedDataFrame(newphenodata,
    varMetadata=varMetadata(object))

  if (inherits(object, "QDNAseqReadCounts")) {
    storage.mode(newcounts) <- "integer"
    object2 <- new("QDNAseqReadCounts", bins=bins, counts=newcounts,
      phenodata=newphenodata)
  } else if (inherits(object, "QDNAseqCopyNumbers")) {
    object2 <- new("QDNAseqCopyNumbers", bins=bins, copynumber=newcopynumber,
      phenodata=newphenodata)
  }
  object2
})

# EOF

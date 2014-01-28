#########################################################################/**
# @RdocFunction poolRuns
#
# @alias poolRuns,QDNAseqReadCounts,character-method
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
setMethod("poolRuns", signature=c(object="QDNAseqReadCounts",
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
  phenodata <- phenoData(object)
  phenodata[, 1] <- samples
  newphenodata <- phenodata[0, ]
  if (class(object) == "QDNAseqReadCounts") {
    counts <- assayDataElement(object, "counts")
    newcounts <- matrix(NA_integer_, nrow=nrow(counts),
      ncol=length(newsamples), dimnames=list(rownames(counts), newsamples))
  } else if (class(object) == "QDNAseqCopyNumbers") {
    copynumber <- assayDataElement(object, "copynumber")
    newcopynumber <- matrix(NA_real_, nrow=nrow(counts),
      ncol=length(newsamples), dimnames=list(rownames(counts), newsamples))
  }
  
  concatenateIfNotEqual <- function(x) {
    x <- sort(unique(x))
    paste(x, collapse=";")
  }

  for (newsample in newsamples) {
    replicates <- samples == newsample
    if (class(object) == "QDNAseqReadCounts") {
      newcounts[, newsample] <- rowSums(counts[, replicates, drop=FALSE])
    } else if (class(object) == "QDNAseqCopyNumbers") {
      newcopynumber[, newsample] <- rowMeans(copynumber[, replicates,
        drop=FALSE])
    }
    oldphenodata <- phenodata[replicates, ]
    totalReads <- sum(oldphenodata$reads)
    numericCols <- sapply(oldphenodata, is.numeric)
    oldphenodata[1, numericCols] <- apply(oldphenodata[, numericCols,
      drop=FALSE], 2, mean)
    oldphenodata[1, !numericCols] <- apply(oldphenodata[, !numericCols,
      drop=FALSE], 2, concatenateIfNotEqual)
    oldphenodata[1, "reads"] <- totalReads
    newphenodata <- rbind(newphenodata, oldphenodata[1,])
  }
  rownames(newphenodata) <- newphenodata[, 1]

  if (class(object) == "QDNAseqReadCounts") {
    object2 <- new("QDNAseqReadCounts", bins=bins, counts=newcounts,
      phenodata=newphenodata)
  } else if (class(object) == "QDNAseqCopyNumbers") {
    object2 <- new("QDNAseqCopyNumbers", bins=bins, copynumber=newcopynumber,
      phenodata=newphenodata)
  }
  object2
})

# EOF

segmentDataWithWeights <- function (input, weights=TRUE, tgr=NULL, ...) 
{
  if (length(weights)==1 & weights & !is.null(tgr)) {
    input <- input[!is.na(tgr),]
    tgr <- tgr[!is.na(tgr)]
    tgr <- abs(tgr)
    tgr[tgr==0] <- min(tgr[tgr!=0], na.rm=TRUE)
    weights <- 1/tgr
  }
  CNA.object <- DNAcopy::CNA(copynumber(input), chromosomes(input), bpstart(input), data.type="logratio")
  cat("Start data segmentation .. \n")
  if (length(weights)==1 && !weights) {
    segmented <- segment(CNA.object, ...)
  } else {
    segmented <- segment(CNA.object, weights=weights, ...)
  }
  numclone <- segmented$output$num.mark
  smrat <- segmented$output$seg
  numsmrat <- cbind(smrat, numclone)
  repdata <- function(row) {
    rep(row[1], row[2])
  }
  makelist <- apply(numsmrat, 1, repdata)
  joined <- unlist(makelist)
  rm(makelist)
  joined <- matrix(joined, ncol=ncol(input), byrow=FALSE)
  joined <- CGHcall:::.assignNames(joined, input)
  result <- CGHcall:::.segFromRaw(input, joined)
  pData(result) <- pData(input)
  result
}

# EOF

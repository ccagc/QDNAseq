#########################################################################/**
# @RdocFunction clonality.analysis 
#
# @alias clonality.analysis,QDNAseqCopyNumbers-method
#
# @title "Clonality testing using copy number data"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \value{
#     Returns an object of class QDNAseqCopyNumbers with segmentation results
#         added.
# }
#
# \examples{
# 
# }
# \references{
#Ostrovnaya, I., Olshen, A. B., Seshan, V.E., Orlow, I., Albertson, D. G. and Begg, C. B. (2010), A metastasis or a second independent cancer? Evaluating the clonal origin of tumors using array copy number data. Statistics in Medicine, 29: 1608-1621
#
#Ostrovnaya, I. and Begg, C.  Testing Clonal Relatedness of Tumors Using Array Comparative Genomic Hybridization: A Statistical Challenge Clin Cancer Res March 1, 2010 16:1358-1367
#
#Venkatraman, E. S. and Olshen, A. B. (2007). A faster circular binary segmentation algorithm for the analysis of array CGH data. Bioinformatics, 23:657 63.
#
#Olshen, A. B., Venkatraman, E. S., Lucito, R., Wigler, M. (2004). Circular binary segmentation for the analysis of array-based DNA copy number data. Biostatistics 5: 557-572.
#
#}
# @author "DS"
#
# \seealso{
#     Internally, @see "Clonality::clonality.analysis" of the \pkg{Clonality} package is used.
# }
#
# @keyword clonality
#*/#########################################################################
#setMethod("likelihoodRatio", signature=c(object="QDNAseqCopyNumbers"),
#          definition=function(obj, patients, maxRow=15000, refdata=NULL, ...) {

likelihoodRatio <- function(obj, patients, maxRow=15000, refdata=NULL, ...) {
            # Prepare data
            
            fData(obj) -> fd
            pData(obj) -> pd
            
            cntrmrs <- CGHbase:::.getCentromere(37)
            chrn <- fd$chromosome
            chrn[chrn == "X"] <- 23
            chrn[chrn == "Y"] <- 24
            
            fd$start / cntrmrs[ as.integer(chrn) ] -> centFact
            
            arm <- rep("q", length(centFact))
            arm[centFact < 1] <- "p"
            
            data <- assayDataElement(obj, "copynumber")
            data <- log2(data[fd$use,])
            
            chrom <- paste("chr", sprintf("%02d", as.integer(chrn[fd$use])), arm[fd$use], sep="")
            maploc <- fd$start[fd$use] / 1000 
            CNA <- CNA(genomdat=data, chrom=chrom, maploc=maploc, data.type="logratio")
            
            if (!is.null(refdata)) {
              reflog2 <- log2(assayDataElement(refdata, "copynumber"))
              reflog2 <- reflog2[fd$use,]
              refdata <- CNA(genomdat=reflog2, chrom=chrom, maploc=maploc, data.type="logratio")
            }
            
            if (!is.na(maxRow) & nrow(data) > maxRow) {
              print("Number of probes exceeds maxRow, averaging")
              k <- (nrow(data) %/% maxRow) + 1
              ave.adj.probes(CNA, k) -> CNA
              if (!is.null(refdata)) {
                ave.adj.probes(refdata, k) -> refdata
              }
            }
            
            clonality.analysis(CNA, patients, refdata=refdata, ...)
}

correlation <- function(obj, labels=NULL, groups=NULL,...) {
  pd <- pData(obj)
  fd <- fData(obj)
  
  if (is.null(labels)) {
    labels <- pd$name
  }
  
  if (!"segmented" %in% assayDataElementNames(obj)) {
    vmsg("Segmented data is required, but not available,... exiting")
    return()
  }
  
  sn <- assayDataElement(obj, "segmented")
  sn <- sn[fd$use,]
  dm <- as.dist(1 - cor(sn))
  hc <- hclust(dm)
  hc$labels <- labels
  hc$groups <- groups
  
  # Distance matric and cluster object
  list(sn=sn, dm=dm, hc=hc)
}

clonality.test <- function(obj = NULL, corData = NULL, llrData = NULL, sbjctLst=NULL, labels=NULL, refdata=NULL, ...) {
  # Do correlation
  if (is.null(corData)) {
   corData <- correlation(obj, labels=labels)
  }
  
  # Do likelihoodRatio
  if (is.null(llrData)) {
    llrData <- likelihoodRatio(obj, patients=sbjctLst, refdata=refdata)
  }
  
  # Combine table
  patients <- sbjctLst
  
  # All combinations
  allPairs <- apply(combn(1:length(patients), 2), 2, function(x) { paste(sprintf("%03d",x), collapse=":")  })
  
  sapply(unique(patients), function(x) {
    cb <- combn(which(patients == x), 2)
    cbStr <- apply(cb, 2, function(x) {paste(sprintf("%03d", x), collapse=":")} )
  }) -> pairs
  
  # clonTab
  pat <- rep(unique(patients), sapply(pairs, length))
  res <- unlist(sapply(pairs, function(x) { corData$dm[ allPairs %in% x  ] }), use.names=F)

  clonTab <- data.frame(patient = pat, combination = unlist(pairs, use.names=F), cor = res, llr2 = log(llrData$LR$LR2), llr2.p = llrData$LR$LR2pvalue, stringsAsFactors=F)
  
  # refTab
  # Correlation
  refCor <- corData$dm[ !allPairs %in% unlist(pairs) ]
  refCorName <- allPairs[ !allPairs %in% unlist(pairs) ]
  
  # LR2
  # Reorder LLR2 data
  refLlr2Name <- paste(llrData$refLR$Sample1, llrData$refLR$Sample2, sep=":")
  refLlr2Name2 <- unlist(sapply(strsplit(refLlr2Name, "\\.|:"), function(x) { paste(sprintf("%03d", sort(as.numeric(x[c(2,4)]))), collapse=":")  } ))
  
  refLR2ord <- order(unlist(refLlr2Name2))
  refLlr2 <- log(llrData$refLR$LR2)
  
  refTab <- data.frame(name = refCorName, check = refLlr2Name[refLR2ord], cor = refCor, llr2 = refLlr2[refLR2ord], stringsAsFactors=F)
  
  # Make sure data types are correct:
  
  # Return data
  # TODO insert into object?
  list(corData = corData, llrData = llrData, clonTab = clonTab, refTab = refTab)
}

# Plot correlation cluster
plot.cor.cluster <- function(corData, scale=T, main="Cluster plot", labwidth=.3, ...) {
  sn <- corData$sn
  hc <- corData$hc
  
  if (scale) {
    sn <- scale(sn)
  }
  
  # Make colors, 50 shades of blue to black to red
  col <- c(rgb(0, 0, seq(1,0, len=25)), rgb(seq(0,1, len=25), 0,0))
  
  # Determine chrom lengths
  chroms <- sapply(strsplit(rownames(sn), ":"), function(x) x[1])
  rle(chroms) -> cc
  
  # Plot
  lmat <- matrix(c(2,3,1), ncol=3)
  lwidth <- c(1,labwidth,5)
  layout(lmat, widths = lwidth)
  par(mar=c(5,0,4,2) + .1)
  snO <- sn[,hc$order]
  image(snO, col=col, zlim=c(-3,3), xaxt="n", yaxt="n", main=main, xlab="chromosome", ...)
  image(snO > 3, col=6, zlim=c(1,1), add=T)
  image(snO < -3, col=5, zlim=c(1,1), add=T)
  
  abline(v=cumsum(cc$lengths) / nrow(sn), lty=2, lwd=.5, col="lightgrey")
  axis(1, at=(cumsum(cc$lengths) - cc$lengths / 2) / nrow(sn), labels=1:22)
  par(mar=c(5, 2, 4, -0.1)  + .1)
  plot(as.dendrogram(hc), horiz=T, yaxs="i", axes=F, leaflab="n")
  
  par(mar=c(5,-0.1,4,-0.1) + .1)
  labels <- hc$labels
  plot(NA, xlim=c(0,1), ylim=c(0, 1), type="n", yaxs="i", axes=F, xlab="patient")
  nl <- length(labels)  
  rect(0, 1:nl / nl, 1, (1:nl -1) / nl, border="white", col="grey")
  text(0, (1:nl - .6) / nl, labels[hc$order], pos=4) 
  par(mar=c(5,4,4,2) + .1)

}

clonality.report <- function(clone.test = NULL, clonTab = NULL, correlations = NULL, LLR2 = NULL, 
  labels = NULL, corThresh=0.54, llr2Thresh=0, refTab = NULL, refLLR2 = NULL, refCorrelations = NULL, main=NULL) {
  
  if (!is.null(clone.test)) {
    clonTab <- clone.test$clonTab
    refTab <- clone.test$refTab
  }
  
  if (!is.null(clonTab)) {
    LLR2 <- clonTab$llr2
    correlations <- clonTab$cor
    labels <- clonTab$patient
  }

  if (!is.null(refTab)) {
    refLLR2 <- refTab$llr2
    refCorrelations <- refTab$cor
  }
  
  if (is.null(labels)) {
    labels <- 1:length(LLR2)
  }
  
  LLR2 >= llr2Thresh -> LRsel
  1 - correlations >= corThresh -> corSel
  
  plot(x="", y="", xlim=c(-40,100), ylim=c(-1,1), ylab="Correlation of segments", xlab="Log likelihood ratio", main=main, bg="lightyellow")
  rect(xleft=c(-50, -5, llr2Thresh), ybottom=c(-2, .45, corThresh), xright=110, ytop=1.3,col=c("lightyellow", "grey", "lightblue"), border=NA)
  abline(v=llr2Thresh, h=corThresh, lty=2)
    
  pch <- rep(16, length(LLR2))
  pch[LRsel & !corSel] <- 21  
  pch[!LRsel & corSel] <- 21
  pch[!LRsel & !corSel] <- 4

  legendLab <- paste(1:length(labels), labels, sep=": ")
  legendPch <- pch
  legendCol <- 1
  
  if (!is.null(refLLR2) & !is.null(refCorrelations)) {
    points(x=refLLR2, y=1 - refCorrelations, pch=4, col="#ff000030", cex=.5)
    legendLab <- c(legendLab, "ref")
    legendPch <- c(legendPch, 4)
    legendCol <- c(rep(1, length(labels)), 2)
  }
  
  points(x=LLR2, y=1 - correlations, pch=pch)
  text(x=LLR2, y=1 - correlations, pos=4, labels=1:length(labels), cex=.8)
  
  legend("bottomright", legend=legendLab, bg="white", pch=legendPch, cex=0.8, col=legendCol)   
}

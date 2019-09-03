#########################################################################/**
# @RdocFunction createBins
#
# @alias calculateMappability
# @alias calculateBlacklist
# @alias calculateBlacklistByRegions
# @alias iterateResiduals
#
# @title "Builds bin annotation data for a particular bin size"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{bsgenome}{A BSgenome package.}
#     \item{binSize}{A @numeric scalar specifying the width of the bins
#         in units of kbp (1000 base pairs), e.g. \code{binSize=15} corresponds
#         to 15 kbp bins.}
#     \item{ignoreMitochondria}{Whether to ignore the mitochondria, defined as
#          chromosomes named 'chrM', 'chrMT', 'M', or 'MT'.}
#     \item{excludeSeqnames}{Character vector of seqnames which should be
#         ignored.}
# }
#
# \value{
#     Returns a @data.frame with columns \code{chromosome, start, end, bases,
#     and gc}, which correspond to the chromosome name, positions of the first
#     and last base pair in the bin, the percentage of characterized nucleotides
#     (A, C, G, or T, i.e. non-N), and GC content (percentage of C and G
#     nucleotides of non-N nucleotides).
# }
#
# \section{Parallel processing}{
# The \pkg{future} is used parallelize the following functions:
#  \itemize{
#   \item \code{createBins()} - parallelizes binned GC content across chromosomes
#   \item \code{calculateBlacklist()} - parallelizes overlap counts across bins)
#  }
# }
#
# \examples{
# \dontrun{# NOTE: These take a very long time to run.
# library(BSgenome.Hsapiens.UCSC.hg19)
# bins <- createBins(BSgenome.Hsapiens.UCSC.hg19, 15)
# bins$mappability <- calculateMappability(bins,
#     bigWigFile='/path/to/wgEncodeCrgMapabilityAlign50mer.bigWig',
#     bigWigAverageOverBed='/path/to/bigWigAverageOverBed')
# bins$blacklist <- calculateBlacklist(bins,
#     bedFiles=c('/path/to/wgEncodeDacMapabilityConsensusExcludable.bed',
#     '/path/to/wgEncodeDukeMapabilityRegionsExcludable.bed'))
# bins$residual <- iterateResiduals(readCountsG1K)
# }
# }
#
# @author "IS"
#
# \seealso{
#     @see "getBinAnnotations".
# }
#*/#########################################################################
createBins <- function(bsgenome, binSize, ignoreMitochondria=TRUE,
    excludeSeqnames=NULL) {
    chrs <- GenomeInfoDb::seqnames(bsgenome)
    try({
        info <- GenomeInfoDb::genomeStyles(GenomeInfoDb::organism(bsgenome))
        style <- GenomeInfoDb::seqlevelsStyle(bsgenome)[1]
        chrs <- info[, style]
    }, silent=TRUE)
    if (!is.null(excludeSeqnames)) {
        chrs <- chrs[!chrs %in% excludeSeqnames]
    }
    if (ignoreMitochondria) {
        selectedMT <- grep("^(chr)?M(T)?$", chrs)
        if (length(selectedMT) != 0L)
            chrs <- chrs[-selectedMT]
    }
    lengths <- GenomeInfoDb::seqlengths(bsgenome)[chrs]
    vmsg("Creating bins of ", binSize, " kbp for genome ",
        substitute(bsgenome))

    # Bin size in units of base pairs
    binWidth <- as.integer(binSize * 1000L)

    chrData <- flapply(chrs, function(chr) {
        vmsg("    Processing ", chr, " ...")
        chr.size <- lengths[chr]
        chr.starts <- seq(from=1L, to=chr.size, by=binWidth)
        chr.ends <- chr.starts + binWidth - 1L
        chr.ends[length(chr.ends)] <- chr.size
        chr.seq <- BSgenome::getSeq(bsgenome, chr, as.character=TRUE)
        bin.seq <- substring(chr.seq, first=chr.starts, last=chr.ends)
        acgt <- gsub("[^ACGT]", "", bin.seq)
        cg <- gsub("[^CG]", "", acgt)
        chr.bases <- nchar(acgt) / (binWidth) * 100
        chr.gc <- nchar(cg) / nchar(acgt) * 100
        list(start=chr.starts, end=chr.ends, bases=chr.bases, gc=chr.gc)
    })
    start <- unlist(lapply(chrData, FUN="[[", "start"))
    end <- unlist(lapply(chrData, FUN="[[", "end"))
    bases <- unlist(lapply(chrData, FUN="[[", "bases"))
    gc <- unlist(lapply(chrData, FUN="[[", "gc"))
    gc[is.nan(gc)] <- NA_real_
    bins <- data.frame(chromosome=rep(chrs, times=ceiling(lengths/binWidth)),
        start, end, bases, gc, stringsAsFactors=FALSE)
    bins$chromosome <- sub("^chr", "", bins$chromosome)
    rownames(bins) <- sprintf("%s:%i-%i", bins$chromosome, bins$start, bins$end)
    bins
}

calculateMappability <- function(bins, bigWigFile,
    bigWigAverageOverBed="bigWigAverageOverBed", chrPrefix="chr") {
    
    binbed <- tempfile(fileext=".bed")
    mapbed <- tempfile(fileext=".bed")
    bins <- bins[, c("chromosome", "start", "end")]
    if (!is.null(chrPrefix) && !is.na(chrPrefix[1]) && chrPrefix[1] != "") {
        vmsg("Prefixing chromosome names with '", chrPrefix[1], "'.")
        bins$chromosome <- paste0(chrPrefix[1], bins$chromosome)
    }
    vmsg("Calculating mappabilities per bin for chromosomes:\n    ",
        paste(unique(bins$chromosome), collapse=", "),
        "\nfrom file:\n    ", bigWigFile, "\nchromosomes to process:   ",
        rep(".", times=length(unique(bins$chromosome))), "\n    ", appendLF=FALSE)
    bins$start <- bins$start - 1
    bins$name <- seq_len(nrow(bins))
    scipen <- options("scipen")
    options(scipen=10)
    write.table(bins, binbed, quote=FALSE, sep="\t", row.names=FALSE,
        col.names=FALSE)
    options(scipen=scipen)
    cmd <- paste0(bigWigAverageOverBed, ' "', bigWigFile, '" "', binbed,
        '" -bedOut="', mapbed, '" /dev/null')
    system(cmd)
    map <- read.table(mapbed, sep="\t", as.is=TRUE)
    map <- map[order(map$V4), ]
    map$V5 * 100
}

calculateBlacklist <- function(bins, bedFiles, ...) {
    vmsg("Calculating overlaps per bin with BED files \n    ", paste(bedFiles,
        collapse="\n    "), "\n    ...", appendLF=FALSE)

    ## Detect deprecated usage of argument 'ncpus'
    args <- list(...)
    if ("ncpus" %in% names(args)) {
      .Deprecated(msg=paste("Argument 'ncpus' of calculateBlacklist() is",
          "deprecated and ignored. Use options(mc.cores=ncpu) instead."))
    }

    beds <- list()
    for (bed in bedFiles)
        beds[[bed]] <- read.table(bed, sep="\t", as.is=TRUE)
    combined <- beds[[1L]]
    if (length(beds) >= 2L)
        for (i in 2:length(beds))
            combined <- rbind(combined, beds[[i]])
    combined <- combined[, 1:3]
    colnames(combined) <- c("chromosome", "start", "end")
    combined$chromosome <- sub("^chr", "", combined$chromosome)
    combined <- combined[combined$chromosome %in% unique(bins$chromosome), ]
    combined <- combined[!is.na(combined$chromosome), ]
    combined$start <- combined$start + 1
    ## define correct sorting order of chromosomes as the order in which they
    ## are in the bins
    chromosomes <- unique(bins$chromosome)
    chromosomeOrder <- factor(combined$chromosome, levels=chromosomes,
       ordered=TRUE)
    combined <- combined[order(chromosomeOrder, combined$start), ]
    joined <- data.frame()
    prev <- combined[1L,]
    # Sanity check
    stopifnot(nrow(combined) >= 2L);
    for (i in 2:nrow(combined)) {
        if (combined[i, "chromosome"] != prev$chromosome ||
            combined[i, "start"] > (prev$end + 1)) {
            joined <- rbind(joined, prev)
            prev <- combined[i,]
        } else {
            prev$end <- max(prev$end, combined[i, "end"])
        }
    }
    joined <- rbind(joined, prev)
    overlap.counter <- function(x, joined) {
        chr <- x["chromosome"]
        start <- as.integer(x["start"])
        end <- as.integer(x["end"])
        overlaps <- joined[joined$chromosome == chr &
                           joined$start      <= end &
                           joined$end        >= start, ]
        bases <- 0
        for (i in rownames(overlaps))
            bases <- bases + min(end, overlaps[i, "end"]) -
                max(start, overlaps[i, "start"]) + 1
        bases / (end - start + 1) * 100
    }
    blacklist <- fapply(bins, MARGIN=1L, FUN=overlap.counter, joined)
    vmsg()
    blacklist
}

iterateResiduals <- function(object, adjustIncompletes=TRUE,
    cutoff=4.0, maxIter=30, ...) {
    first <- sum(binsToUse(object))
    previous <- first
    vmsg("Iteration #1 with ", format(previous, big.mark=","),
        " bins.")
    object <- estimateCorrection(object, ...)
    counts <- assayDataElement(object, "counts")
    if (adjustIncompletes) {
        counts <- counts / fData(object)$bases * 100L
        counts[fData(object)$bases == 0] <- 0L
    }
    fit <- assayDataElement(object, "fit")
    residuals <- counts / fit - 1
    # residuals[fit == 0] <- NA_real_
    residuals[!binsToUse(object), ] <- NA_real_
    residual <- rowMedians(residuals, na.rm=TRUE)
    cutoffValue <- cutoff * madDiff(residual, na.rm=TRUE)
    if (is.numeric(cutoff))
        binsToUse(object) <- binsToUse(object) & !is.na(residual) &
            abs(residual) <= cutoffValue
    num <- sum(binsToUse(object))
    iter <- 2
    while (previous != num && iter <= maxIter) {
        previous <- num
        vmsg("Iteration #", iter, " with ", format(previous, big.mark=","),
            " bins.")
        object <- estimateCorrection(object, ...)
        fit <- assayDataElement(object, "fit")
        residuals <- counts / fit - 1
        # residuals[fit == 0] <- NA_real_
        residuals[!binsToUse(object), ] <- NA_real_
        residual <- rowMedians(residuals, na.rm=TRUE)
        binsToUse(object) <- binsToUse(object) & !is.na(residual) &
            abs(residual) <= cutoffValue
        num <- sum(binsToUse(object))
        iter <- iter + 1
    }
    if (previous == num) {
        vmsg("Convergence at ", format(previous, big.mark=","), " bins.")
        vmsg(format(first-previous, big.mark=","), " additional bins removed.")
    } else if (iter == maxIter) {
        vmsg("Reached maxIter=", maxIter, " iterations without convergence.")
    }
    residual
}


###############################################################################
# 
# Calculates percentage overlap with regions <- cbind(chromosome, start, end)
# More generic, eg when bed files are not readily available.
# 
###############################################################################
calculateBlacklistByRegions <- function(bins, regions) {
    vmsg("Calculating percent overlap per bin with regions")
    
    combined <- as.data.frame(regions)
    
    colnames(combined) <- c("chromosome", "start", "end")
    # Remove chr if present
    combined$chromosome <- sub("^chr", "", combined$chromosome)
    
    # Remove chromosome not pressent in bin
    combined <- combined[combined$chromosome %in% unique(bins$chromosome), ]
    
    # Convert XY to number integer
    combined$chromosome[combined$chromosome=="X"] <- "23"
    combined$chromosome[combined$chromosome=="Y"] <- "24"
    combined$chromosome <- as.integer(combined$chromosome)
    combined <- combined[!is.na(combined$chromosome), ]
    
    bins$chromosome[bins$chromosome=="X"] <- "23"
    bins$chromosome[bins$chromosome=="Y"] <- "24"
    bins$chromosome <- as.integer(bins$chromosome)
    
    # Assume 1 based coordinate system
    # TODO Check wether necessary
    
    tmp <- combined$start + 1
    combined$start <- tmp
    combined <- combined[order(combined$chromosome, combined$start), ]
    
    # Determine binsize
    binSize <- diff(bins$start[1:2])
    
    # Calculate index, residual, 
    # add 1 because it serves as idx
    combined$si <- combined$start %/% binSize + 1 
    combined$sm <- binSize - combined$start %% binSize
    # add 1 because it serves as idx
    combined$ei <- combined$end %/% binSize + 1 
    combined$em <- combined$end %% binSize
    combined$seDiff <- combined$end - combined$start 
    combined$idDiff <- combined$ei - combined$si 
    
    # Calculate continuous IDX
    c(0, cumsum(rle(bins$chromosome)$lengths)) -> chrI
    combined$sI <- combined$si + chrI[ as.integer(combined$chromosome) ]
    combined$eI <- combined$ei + chrI[ as.integer(combined$chromosome) ]
    
    vmsg("Processing partial overlaps")
    # Partial overlaps of segments
    sel1 <- combined$idDiff >= 1
    # Sum complete overlaps eg segment smaller than bin
    sel2 <- combined$idDiff == 0
    
    aggregate(c(
        combined$sm[sel1], 
        combined$em[sel1],
        combined$seDiff[sel2]
    ), 
        list(c(
            combined$sI[sel1], 
            combined$eI[sel1],
            combined$sI[sel2]
        )
        ), max) -> res12
    
    vmsg("Complete overlaps")
    # Sum complete overlaps of segments eg segment larger than bin
    sel3 <- combined$idDiff > 1
    unlist(sapply(which(sel3), FUN=function(x) {
        s <- combined$sI[x] + 1
        e <- combined$eI[x] - 1
        s:e
    })) -> res3
    
    res <- rbind(res12, cbind(Group.1 = res3, x = rep(binSize, times=length(res3))))
    
    aggregate(res$x, list(res$Group.1), max) -> res
    
    res$x / binSize * 100 -> res$pct
    
    blacklist <- rep(0, times=nrow(bins))
    blacklist[ as.numeric(res$Group.1) ] <- res$pct
    
    blacklist
}



# EOF

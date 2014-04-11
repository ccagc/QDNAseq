#########################################################################/**
# @RdocFunction createBins
#
# @alias calculateMappability
# @alias calculateBlacklist
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
#     \item{ignoreMitochondria}{Wheter to ignore the mitochondria.}
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
#     @see "downloadBinAnnotations".
# }
#*/#########################################################################
createBins <- function(bsgenome, binSize, ignoreMitochondria=TRUE) {
    info <- GenomeInfoDb::genomeStyles(organism(bsgenome))
    style <- GenomeInfoDb::seqlevelsStyle(bsgenome)
    chrs <- info[, style]
    if (ignoreMitochondria)
        chrs <- chrs[-grep("^(chr)?M(T)?$", chrs)]
    lengths <- GenomicRanges::seqlengths(bsgenome)[chrs]
    start <- end <- integer()
    bases <- gc <- numeric()
    vmsg("Creating bins of ", binSize, " kbp for genome ",
        substitute(bsgenome))

    # Bin size in units of base pairs
    binWidth <- binSize*1000L

    for (chr in chrs) {
        vmsg("    Processing ", chr, " ...", appendLF=FALSE)
        chr.size <- lengths[chr]
        chr.starts <- seq(from=1, to=chr.size, by=binWidth)
        chr.ends <- chr.starts + binWidth - 1L
        chr.ends[length(chr.ends)] <- chr.size
        chr.seq <- BSgenome::getSeq(bsgenome, chr, as.character=TRUE)
        bin.seq <- substring(chr.seq, first=chr.starts, last=chr.ends)
        acgt <- gsub("[^ACGT]", "", bin.seq)
        cg <- gsub("[^CG]", "", acgt)
        chr.bases <- nchar(acgt) / (binWidth) * 100
        chr.gc <- nchar(cg) / nchar(acgt) * 100
        start <- c(start, chr.starts)
        end <- c(end, chr.ends)
        bases <- c(bases, chr.bases)
        gc <- c(gc, chr.gc)
        vmsg()
    }
    gc[is.nan(gc)] <- NA_real_
    bins <- data.frame(chromosome=rep(chrs, times=ceiling(lengths/binWidth)),
        start, end, bases, gc, stringsAsFactors=FALSE)
    bins$chromosome <- sub("^chr", "", bins$chromosome)
    rownames(bins) <- sprintf("%s:%i-%i", bins$chromosome, bins$start, bins$end)
    bins
}

calculateMappability <- function(bins, bigWigFile,
    bigWigAverageOverBed="bigWigAverageOverBed", chrPrefix="chr") {
    vmsg("Calculating mappabilities per bin from file\n    ", bigWigFile,
        "\n    ", appendLF=FALSE)
    binbed <- tempfile(fileext=".bed")
    mapbed <- tempfile(fileext=".bed")
    bins <- bins[,c("chromosome", "start", "end")]
    bins$chromosome <- paste0(chrPrefix, bins$chromosome)
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
    map$V5 * 100
}

calculateBlacklist <- function(bins, bedFiles, ncpus=1) {
    vmsg("Calculating overlaps per bin with BED files \n    ", paste(bedFiles,
        collapse="\n    "), "\n    ...", appendLF=FALSE)
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
    combined$chromosome[combined$chromosome=="X"] <- "23"
    combined$chromosome[combined$chromosome=="Y"] <- "24"
    combined$chromosome <- as.integer(combined$chromosome)
    combined <- combined[!is.na(combined$chromosome), ]
    combined$start <- combined$start + 1
    combined <- combined[order(combined$chromosome, combined$start), ]
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
    bins$chromosome[bins$chromosome=="X"] <- "23"
    bins$chromosome[bins$chromosome=="Y"] <- "24"
    bins$chromosome <- as.integer(bins$chromosome)
    overlap.counter <- function(x, joined) {
        chr <- as.integer(x["chromosome"])
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
    if (ncpus > 1L) {
        snowfall::sfInit(parallel=TRUE, cpus=ncpus)
        snowfall::sfExport(list=c("overlap.counter"))
        blacklist <- snowfall::sfApply(bins, 1, overlap.counter, joined)
        snowfall::sfStop()
    } else {
        blacklist <- apply(bins, MARGIN=1L, FUN=overlap.counter, joined)
    }
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
    # residuals[fit == 0] <- NA
    residuals[!binsToUse(object), ] <- NA
    residual <- apply(residuals, 1, median, na.rm=TRUE)
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
        # residuals[fit == 0] <- NA
        residuals[!binsToUse(object), ] <- NA
        residual <- apply(residuals, 1, median, na.rm=TRUE)
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

# EOF

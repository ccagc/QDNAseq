#########################################################################/**
# @RdocFunction exportBins
#
# @alias exportBins,QDNAseqSignals-method
#
# @title "Exports to a file"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{object}{A @see "QDNAseqReadCounts" or @see "QDNAseqCopyNumbers"
#         object.}
#     \item{file}{Filename. For formats that support only one sample per file,
#         such as BED, '\%s' can be used as a placeholder for sample name or
#         '\%d' for sample number.}
#     \item{format}{Format to export in. Currently supported ones are "tsv" (tab
#         separated values), "igv" (Integrative Genomics Viewer), and "bed" (BED
#         file format).}
#     \item{type}{Type of data to export, options are "copynumber" (corrected or
#         uncorrected read counts), "segments", or "calls".}
#     \item{filter}{If @TRUE, bins are filtered, otherwise not.}
#     \item{logTransform}{If @TRUE (default), data will be log2-transformed.}
#     \item{digits}{The number of digits to round to. If not @numeric, no
#         no rounding is performed.}
#     \item{chromosomeReplacements}{A named character vector of chromosome name
#         replacements to be done. Only used when \code{object} is of class
#         @see "cghRaw", @see "cghSeg", @see "cghCall", or @see "cghRegions",
#        since these classes store chromosome names as integers, whereas all
#        QDNAseq object types use character vectors. Defaults to
#         \code{c("23"="X", "24"="Y", "25"="MT")} for human.}
#     \item{...}{Additional arguments passed to @see "utils::write.table".}
# }
#
# \details{
#     Exports \code{object} to a file.
# }
#
# \examples{
# \dontrun{
# data(LGG150)
# readCounts <- LGG150
# readCountsFiltered <- applyFilters(readCounts)
# readCountsFiltered <- estimateCorrection(readCountsFiltered)
# copyNumbers <- correctBins(readCountsFiltered)
# copyNumbersNormalized <- normalizeBins(copyNumbers)
# copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
# exportBins(copyNumbersSmooth, file="LGG150.igv", format="igv")
# }
# }
#
# @author "IS"
#
# @keyword IO
# @keyword file
#*/#########################################################################
exportBins <- function(object, file, 
    format=c("tsv", "igv", "bed", "vcf", "seg"),
    type=c("copynumber", "segments", "calls"),
    filter=TRUE, logTransform=TRUE, digits=3,
    chromosomeReplacements=c("23"="X", "24"="Y", "25"="MT"), ...) {

    makeFilenames <- function(file, names) {
        if (length(grep("%s", file) > 0L)) {
            files <- sapply(seq_along(names), FUN = function(idx) {
                sprintf(file, names[idx])
            })
        } else if (length(grep("%i", file) > 0L)) {
            files <- sapply(seq_along(names), FUN = function(idx) {
                sprintf(file, idx)
            })
        } else {
            if (length(names) > 1L) {
              stop("Argument 'file' must be an sprintf-formatted strings with either a %s or a %i component: ", sQuote(file))
            }
            files <- file
        }
        files
    } # makeFilenames()

    format <- match.arg(format)
    type <- match.arg(type)
    
    if (inherits(object, "QDNAseqSignals")) {
        if (filter) {
            object <- object[binsToUse(object), ]
        }
        feature <- featureNames(object)
        chromosome <- fData(object)$chromosome
        start <- fData(object)$start
        end <- fData(object)$end
        if (inherits(object, "QDNAseqReadCounts")) {
            if (type != "copynumber")
                warning("Ignoring argument 'type' and returning read counts.")
            dat <- assayDataElement(object, "counts")
            type <- "read counts"
        } else {
            if (type == "copynumber") {
                dat <- assayDataElement(object, "copynumber")
            } else if (type == "segments") {
                if (!"segmented" %in% assayDataElementNames(object))
                    stop("Segments not found, please run segmentBins() first.")
                dat <- assayDataElement(object, "segmented")
            } else if (type == "calls") {
                if (!"calls" %in% assayDataElementNames(object))
                    stop("Calls not found, please run callBins() first.")
                dat <- assayDataElement(object, "calls")
            }
        }
        if (logTransform) {
            dat <- log2adhoc(dat)
        }
    } else if (inherits(object, c("cghRaw", "cghSeg", "cghCall",
        "cghRegions"))) {

        feature <- featureNames(object)
        chromosome <- as.character(chromosomes(object))
        for (chromosomeReplacement in names(chromosomeReplacements)) {
            chromosome[chromosome == chromosomeReplacement] <-
                chromosomeReplacements[chromosomeReplacement]
        }
        start <- bpstart(object)
        end <- bpend(object)
        if (inherits(object, c("cghRaw", "cghSeg", "cghCall"))) {
            if (type == "copynumber") {
                dat <- copynumber(object)
            } else if (type == "segments") {
                if (!"segmented" %in% assayDataElementNames(object))
                    stop("Segments not found, please run segmentData() first.")
                dat <- segmented(object)
            } else if (type == "calls") {
                if (!"calls" %in% assayDataElementNames(object))
                    stop("Calls not found, please run CGHcall() first.")
                dat <- calls(object)
            }
        } else if (inherits(object, "cghRegions")) {
            if (type != "calls")
                warning("Ignoring argument 'type' and returning calls.")
            dat <- regions(object)
        }
    }
    
    if (is.numeric(digits)) {
        dat <- round(dat, digits=digits)
    }
    
    oopts2 <- options(scipen=15)
    on.exit(options(scipen=oopts2), add=TRUE)
    
    if (format == "tsv") {
        out <- data.frame(feature=feature, chromosome=chromosome, start=start,
            end=end, dat, check.names=FALSE, stringsAsFactors=FALSE)
        write.table(out, file=file,
            quote=FALSE, sep="\t", na="", row.names=FALSE, ...)
        stopifnot(file_test("-f", file))
    } else if (format == "igv") {
        out <- data.frame(chromosome=chromosome, start=start, end=end,
            feature=feature, dat, check.names=FALSE, stringsAsFactors=FALSE)
        cat('#type=COPY_NUMBER\n#track coords=1\n', file=file)
        suppressWarnings(write.table(out, file=file, append=TRUE,
            quote=FALSE, sep="\t", na="", row.names=FALSE, ...))
        stopifnot(file_test("-f", file))
    } else if (format == "bed" ) {
        names <- sampleNames(object)
        files <- makeFilenames(file, names = names)
        for (i in seq_along(names)) {
            out <- data.frame(chromosome=chromosome, start=start-1, end=end,
                feature=feature, dat[, i, drop=FALSE], strand="+",
                check.names=FALSE, stringsAsFactors=FALSE)
            cat('track name="', names[i], '" description="', type, '"\n',
                file=files[i], sep="")
            suppressWarnings({
                write.table(out, file=files[i], append=TRUE, col.names=FALSE,
                            quote=FALSE, sep="\t", na="", row.names=FALSE, ...)
            })
            stopifnot(file_test("-f", files[i]))
        }
    } else if (format == "vcf") {
        names <- sampleNames(object)
        files <- makeFilenames(file, names = names)
	exportVCF(object, fnames = files)
    } else if (format == "seg") {
        names <- sampleNames(object)
        files <- makeFilenames(file, names = names)
        exportSEG(object, fnames = files)
    }
}


exportVCF <- function(obj, fnames) {
    calls <- assayDataElement(obj, "calls")
    segments <- log2adhoc(assayDataElement(obj, "segmented"))

    fd <- fData(obj)
    pd <- pData(obj)

    vcfHeader <- cbind(c(
			 '##fileformat=VCFv4.2',
			 paste('##source=QDNAseq-', packageVersion("QDNAseq"), sep=""),
			 '##REF=<ID=DIP,Description="CNV call">',
			 '##ALT=<ID=DEL,Description="Deletion">',
			 '##ALT=<ID=DUP,Description="Duplication">',
			 '##FILTER=<ID=LOWQ,Description="Filtered due to call in low quality region">',
			 '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of variant: DEL,DUP,INS">',
			 '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of variant">',
			 '##INFO=<ID=BINS,Number=1,Type=Integer,Description="Number of bins in call">',
			 '##INFO=<ID=SCORE,Number=1,Type=Integer,Description="Score of calling algorithm">',
			 '##INFO=<ID=LOG2CNT,Number=1,Type=Float,Description="Log 2 count">', 
			 '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
			 ))

    oopts2 <- options(scipen=100)
    on.exit(options(scipen=oopts2), add=TRUE)

    for (i in 1:ncol(calls)) {	
	d <- cbind(fd[,1:3], calls[,i], segments[,i])
	sel <- d[,4] != 0 & !is.na(d[,4])

	dsel <- d[sel,]

	rleD <- rle(paste(d[sel,1], d[sel,4], sep=":"))

	endI <- cumsum(rleD$lengths)
	posI <- c(1, endI[-length(endI)] + 1)

	chr <- dsel[posI,1]
	pos <- dsel[posI,2]
	end <- dsel[endI,3]
	score <- dsel[posI,4]
	segVal <- round(dsel[posI,5], digits=2)

	svtype <- rep(NA_character_, times=length(chr)) 
	svlen <- rep(NA_real_, times=length(chr)) 
	gt <- rep(NA_character_, times=length(chr)) 
	bins <- rleD$lengths
	svtype[dsel[posI,4] <= -1] <- "DEL"
	svtype[dsel[posI,4] >= 1] <- "DUP"
	svlen <- end - pos + 1

	gt[score == -2] <- "1/1"	
	gt[score == -1] <- "0/1"	
	gt[score == 1] <- "0/1"	
	gt[score == 2] <- "0/1"	
	gt[score == 3] <- "0/1"	

	id <- "."
	ref <- "<DIP>"
	alt <- paste("<", svtype, ">", sep="")
	qual <- 1000
	filter <- "PASS"
	info <- paste("SVTYPE=", svtype, ";END=", end, ";SVLEN=", svlen, ";BINS=", bins, ";SCORE=", score, ";LOG2CNT=", segVal, sep="")
	format <- "GT"
	sample <- gt
	out <- cbind(chr, pos, id, ref, alt, qual, filter, info, format, sample)	
	colnames(out) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", pd$name[i])

	write.table(vcfHeader, file=fnames[i], quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
	suppressWarnings(write.table(out, file=fnames[i], quote=FALSE, sep="\t", append=TRUE, col.names=TRUE, row.names=FALSE))
        stopifnot(file_test("-f", fnames[i]))
    }
}



exportSEG <- function(obj, fnames=NULL) {

    calls <- assayDataElement(obj, "calls")
    segments <- log2adhoc(assayDataElement(obj, "segmented"))

    fd <- fData(obj)
    pd <- pData(obj)
  
    if (is.null(fnames)) 
	fnames <- pd$name

    if (length(fnames) != length(pd$name)) {
        stop("Length of 'fnames' is too short: ", length(fnames), " != ", length(pd$name))
    }
 
    oopts2 <- options(scipen=100)
    on.exit(options(scipen=oopts2), add=TRUE)

    for (i in 1:ncol(calls)) {	
	d <- cbind(fd[,1:3],calls[,i], segments[,i])
	sel <- d[,4] != 0 & !is.na(d[,4])

	dsel <- d[sel,]

	rleD <- rle(paste(d[sel,1], d[sel,4], sep=":"))

	endI <- cumsum(rleD$lengths)
	posI <- c(1, endI[-length(endI)] + 1)

	chr <- dsel[posI,1]
	pos <- dsel[posI,2]
	end <- dsel[endI,3]
	score <- dsel[posI,4]
	segVal <- round(dsel[posI,5],digits=2)
	bins <- rleD$lengths

	out <- cbind(fnames[i], chr, pos, end, bins, segVal)
	colnames(out) <- c("SAMPLE_NAME", "CHROMOSOME", "START", "STOP", "DATAPOINTS", "LOG2_RATIO_MEAN")

	write.table(out, file = fnames[i], quote=FALSE, sep="\t", append=FALSE, col.names=TRUE, row.names=FALSE)
        stopifnot(file_test("-f", fnames[i]))
    }
}

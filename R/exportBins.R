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
exportBins <- function(object, file, format=c("tsv", "igv", "bed"),
    type=c("copynumber", "segments", "calls"),
    filter=TRUE, logTransform=TRUE, digits=3, ...) {

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
        chromosome[chromosome == "23"] <- "X"
        chromosome[chromosome == "24"] <- "Y"
        chromosome[chromosome == "25"] <- "MT"
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
    tmp <- options("scipen")
    options(scipen=15)
    if (format == "tsv") {
        out <- data.frame(feature=feature, chromosome=chromosome, start=start,
            end=end, dat, check.names=FALSE, stringsAsFactors=FALSE)
        write.table(out, file=file,
            quote=FALSE, sep="\t", na="", row.names=FALSE, ...)
    } else if (format == "igv") {
        out <- data.frame(chromosome=chromosome, start=start, end=end,
            feature=feature, dat, check.names=FALSE, stringsAsFactors=FALSE)
        cat('#type=COPY_NUMBER\n#track coords=1\n', file=file)
        suppressWarnings(write.table(out, file=file, append=TRUE,
            quote=FALSE, sep="\t", na="", row.names=FALSE, ...))
    }  else if (format == "bed" ) {
        for (i in seq_along(sampleNames(object))) {
            if (length(grep("%s", file) > 0L)) {
                filename <- sprintf(file, sampleNames(object)[i])
            } else {
                filename <- sprintf(file, i)
            }
            out <- data.frame(chromosome=chromosome, start=start-1, end=end,
                feature=feature, dat[, i, drop=FALSE], strand="+",
                check.names=FALSE, stringsAsFactors=FALSE)
            cat('track name="', sampleNames(object)[i], '" description="',
                type, '"\n', file=filename, sep="")
            suppressWarnings(write.table(out, file=filename, append=TRUE,
                col.names=FALSE,
                quote=FALSE, sep="\t", na="", row.names=FALSE, ...))
        }
    }
    options(scipen=tmp)
}

# EOF

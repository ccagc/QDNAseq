#########################################################################/**
# @RdocFunction getBinAnnotations
#
# @title "Gets bin annotation data for a particular bin size"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{binSize}{A @numeric scalar specifying the width of the bins
#        in units of kbp (1000 base pairs), e.g. \code{binSize=15} corresponds
#        to 15 kbp bins.}
#     \item{genome}{A @character string specify the genome and genome version
#        to be used.}
#     \item{type}{A @character string specify the experiment type, e.g. "SR50"
#        or "PE100".}
#     \item{path}{A @character string specifying the path for the bin
#         annotation files. Defaults to downloading from the Internet, but can
#         also be a local path. Can also be defined by setting the
#         \code{QDNAseq::binAnnotationPath} option.}
# }
#
# \details{
#     Gets bin annotation data for a particular bin size
# }
#
# \value{
#     Returns a @see "Biobase::AnnotatedDataFrame" object.
# }
#
# \examples{
# \dontrun{
# bins <- getBinAnnotations(15)
# }
# }
#
# @author "IS"
#
# \seealso{
#     @see "createBins".
# }
#
# @keyword IO
#*/#########################################################################
getBinAnnotations <- function(binSize, genome='hg19', type='SR50',
    path=getOption("QDNAseq::binAnnotationPath", NULL)) {

    bins <- NULL

    # when no path argument was provided, try to load from annotation package
    if (is.null(path)) {
        # Check for annotation package
        PKname <- sprintf('QDNAseq.%s', genome)
        PKfound <- PKname %in% .packages(all.available=TRUE)
        if(!PKfound) {
            vmsg("Annotation package ", PKname, " not found.")
            path <- "http://qdnaseq-ireland.s3.amazonaws.com"
        } else {
            # Look for binAnnotations with data()
            BAname <- sprintf('%s.%gkbp.%s', genome, binSize, type)
            tryCatch(data(list=BAname, package=PKname, envir = environment()),
                warning=function(x) {
                    # TODO suggest install package
                    vmsg(BAname, " not found in local datasets")
                }
            )
            try(bins <- get(BAname), silent=TRUE)
            if(!is.null(bins)) {
                vmsg('Loaded bin annotations for genome ', genome,
                    ', bin size ', binSize, 'kbp, and experiment type ', type,
                    ' from annotation package ', PKname, '.')
                return(bins)
            }
        }
    }

    filename <- sprintf('QDNAseq.%s.%gkbp.%s.rds', genome, binSize, type)
    if (length(grep("^http[s]?://", tolower(path[1]))) == 1) {
        vmsg('Downloading bin annotations for genome ', genome,
            ', bin size ', binSize, 'kbp, and experiment type ', type,
            ' from URL ', path, ' ...', appendLF=FALSE)
        remotefile <- file.path(path, filename, fsep='/')
        localfile <- tempfile()
        tryCatch({
            result <- downloadFile(remotefile, localfile)
        }, error=function(e) {
            vmsg(' not found.')
            stop(e)
        })
        bins <- readRDS(localfile)
        file.remove(localfile)
    } else {
        localfile <- file.path(path, filename)
        if (!file.exists(localfile)) {
            stop("File not found: ", localfile)
        }
        vmsg('Loading bin annotations for genome ', genome,
            ', bin size ', binSize, 'kbp, and experiment type ', type,
            ' from path ', path, ' ...', appendLF=FALSE)
        localfile <- file.path(path, filename)
        bins <- readRDS(localfile)
    }
    vmsg()
    bins
}

setMethod("show", signature=c(object="AnnotatedDataFrame"),
    definition=function(object) {
    ## Import private functions
    ns <- asNamespace("Biobase")
    .showAnnotatedDataFrame <- get(".showAnnotatedDataFrame", envir=ns, mode="function")

    if (!is.null(attr(object, "QDNAseqVersion"))) {
        cat("QDNAseq bin annotations\n")
        cat("WARNING: These bin annotations are outdated, please re-download",
        " from URL:\n", "http://cdn.bitbucket.org/ccagc/qdnaseq/downloads/",
        "QDNAseq.hg19.", round((object$end[1] - object$start[1] + 1) / 1000),
        "kbp.SR50.rds\n", sep="")
    } else if (!is.null(attr(object, "QDNAseq"))) {
        info <- attr(object, "QDNAseq")
        cat("QDNAseq bin annotations for ", info$organism,
            ", build ", info$build, ".\n", sep="")
        cat("Created by ", info$author,
            " with QDNAseq ", as.character(info$version),
            ", ", as.character(info$date), ".\n", sep="")
        if (info$version < "0.7.5")
            cat("WARNING: These bin annotations are outdated, ",
                "please re-download from URL:\n",
                "http://cdn.bitbucket.org/ccagc/qdnaseq/downloads/",
                "QDNAseq.hg19.",
                round((object$end[1] - object$start[1] + 1) / 1000),
                "kbp.SR50.rds\n", sep="")
    }
    .showAnnotatedDataFrame((object))
})

# EOF

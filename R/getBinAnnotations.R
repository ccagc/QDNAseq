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
#         annotation file to be downloaded.  The path can either be on the
#         local file system or a URL online.
#         If @NULL (default), then data loaded from an \R package named
#         \pkg{QDNAseq.\{\{genome\}\}}.
#         The default value can be controlled via \R options
#         \code{QDNAseq::binAnnotationPath}.}
#     \item{verbose}{If @TRUE, verbose messages are produced.}
# }
#
# \details{
#     Gets bin annotation data for a particular bin size.
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
getBinAnnotations <- function(binSize, genome="hg19", type="SR50",
    path=getOption("QDNAseq::binAnnotationPath", NULL),
    verbose=getOption("QDNAseq::verbose", TRUE)) {

    oopts <- options("QDNAseq::verbose"=verbose)
    on.exit(options(oopts))

    bins <- NULL

    vmsg_source <- NULL
    
    # Use QDNAseq.{{genome}} annotation package
    if (is.null(path)) {
        # Check for package
        pkgName <- sprintf("QDNAseq.%s", genome)
        pkgPath <- system.file(package = pkgName, mustWork = FALSE)
        if (pkgPath == "") {
          stop("Package not installed: ", pkgName)
        }

        dataPath <- system.file("data", package = pkgName, mustWork = FALSE)
        if (dataPath == "") {
          stop("Package does not provide package data: ", pkgName)
        }

        # Load annotation
        env <- new.env()
        baName <- sprintf("%s.%gkbp.%s", genome, binSize, type)
        tryCatch({
          data(list=baName, package=pkgName, envir=env)
        }, warning = function(w) {
          names <- dir(path=dataPath, pattern="[.]rda$")
          names <- gsub("[.]rda$", "", names)
          stop(sprintf("QDNAseq bin annotation %s was not found in package %s. Available bin annotations are: %s", sQuote(baName), sQuote(pkgName), paste(sQuote(names), collapse=", ")))
        })
        
        bins <- get(baName, mode="S4", envir=env)
        attr(bins, "source") <- sprintf("%s v%s", pkgName, packageVersion(pkgName))
        
        vmsg_source <- paste("annotation package", attr(bins, "source"))
    } else {
      filename <- sprintf("QDNAseq.%s.%gkbp.%s.rds", genome, binSize, type)
  
      ## Download or a local file?
      if (grepl("^http[s]?://", tolower(path[1]))) {
          url <- file.path(path, filename, fsep="/")
          localfile <- tempfile()
          on.exit(file.remove(localfile))
          vmsg("Downloading: ", url)
          result <- downloadFile(url, localfile)
          bins <- readRDS(localfile)
          vmsg_source <- url
      } else {
          localfile <- file.path(path, filename)
          if (!file_test("-f", localfile)) stop("File not found: ", localfile)
          bins <- readRDS(localfile)
          vmsg_source <- localfile
      }
    }
    
    vmsg("Loaded bin annotations for genome ", sQuote(genome), ", bin size ",
         binSize, " kbp, and experiment type ", sQuote(type), " from ",
	 vmsg_source)
 
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

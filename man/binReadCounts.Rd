\name{binReadCounts}
\alias{binReadCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
binReadCounts(binsize, files = NULL, path = ".", ext = "bam", names = NULL, genome = "hg19", cache = TRUE, samtools = "samtools", f = "", F = "0x0404", q = 37, maxChunk = 1e+08)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{binsize}{
%%     ~~Describe \code{binsize} here~~
}
  \item{files}{
%%     ~~Describe \code{files} here~~
}
  \item{path}{
%%     ~~Describe \code{path} here~~
}
  \item{ext}{
%%     ~~Describe \code{ext} here~~
}
  \item{names}{
%%     ~~Describe \code{names} here~~
}
  \item{genome}{
%%     ~~Describe \code{genome} here~~
}
  \item{cache}{
%%     ~~Describe \code{cache} here~~
}
  \item{samtools}{
%%     ~~Describe \code{samtools} here~~
}
  \item{f}{
%%     ~~Describe \code{f} here~~
}
  \item{F}{
%%     ~~Describe \code{F} here~~
}
  \item{q}{
%%     ~~Describe \code{q} here~~
}
  \item{maxChunk}{
%%     ~~Describe \code{maxChunk} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (binsize, files = NULL, path = ".", ext = "bam", names = NULL, 
    genome = "hg19", cache = TRUE, samtools = "samtools", f = "", 
    F = "0x0404", q = 37, maxChunk = 1e+08) 
{
    genome.build <- as.integer(gsub("[^0-9]", "", genome))
    if (genome.build == 36 || genome.build == 18) {
        annfile <- file.path("/scratch/data/ref/hg18/MPScall", 
            paste("QDNAseq.", binsize * 1000, ".txt.gz", sep = ""))
    }
    else {
        annfile <- file.path("/scratch/data/ref/hg19/MPScall", 
            paste("QDNAseq.", binsize * 1000, ".txt.gz", sep = ""))
    }
    if (!file.exists(annfile)) 
        stop("No annotations file found for genome ", genome, 
            " and bin size of ", binsize, "kbp. Please generate it first.")
    bins <- read.table(annfile, header = TRUE, sep = "\t", as.is = TRUE)
    rownames(bins) <- paste(bins$chromosome, ":", bins$start, 
        "-", bins$end, sep = "")
    if (is.null(files)) 
        files <- list.files(path, paste("\\.", ext, "$", sep = ""))
    if (is.null(names)) {
        names <- sub(paste("\\.", ext, "$", sep = ""), "", files)
    }
    else if (length(files) != length(names)) {
        stop("files and names have to be of same length.")
    }
    counts <- matrix(nrow = nrow(bins), ncol = length(names), 
        dimnames = list(rownames(bins), names))
    for (i in 1:length(files)) {
        counts[, i] <- .binReadCountsPerSample(binsize, files[i], 
            path, bins, cache, samtools, f, F, q, maxChunk)
        gc(FALSE)
    }
    dat <- list(bins = bins, counts = counts)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

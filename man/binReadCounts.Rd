\name{binReadCounts}
\alias{binReadCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{binReadCounts
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
binReadCounts(bins, bamfiles = NULL, path = ".", ext = "bam", bamnames = NULL, phenofile = NULL, genome = "hg19", cache = TRUE, samtools = "samtools", f = "", F = "0x0404", q = 37, maxChunk = 100000000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bins}{
%%     ~~Describe \code{bins} here~~
}
  \item{bamfiles}{
%%     ~~Describe \code{bamfiles} here~~
}
  \item{path}{
%%     ~~Describe \code{path} here~~
}
  \item{ext}{
%%     ~~Describe \code{ext} here~~
}
  \item{bamnames}{
%%     ~~Describe \code{bamnames} here~~
}
  \item{phenofile}{
%%     ~~Describe \code{phenofile} here~~
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
function (bins, bamfiles = NULL, path = ".", ext = "bam", bamnames = NULL, 
    phenofile = NULL, genome = "hg19", cache = TRUE, samtools = "samtools", 
    f = "", F = "0x0404", q = 37, maxChunk = 100000000) 
{
    if (is.null(bamfiles)) 
        bamfiles <- list.files(path, paste("\\.", ext, "$", sep = ""))
    if (length(bamfiles) == 0) 
        stop("No files to process.")
    if (is.null(bamnames)) {
        bamnames <- sub(paste("\\.", ext, "$", sep = ""), "", 
            bamfiles)
    }
    else if (length(bamfiles) != length(bamnames)) {
        stop("bamfiles and bamnames have to be of same length.")
    }
    phenodata <- data.frame(name = bamnames, row.names = bamnames, 
        stringsAsFactors = FALSE)
    if (!is.null(phenofile)) {
        pdata <- read.table(phenofile, header = TRUE, sep = "\t", 
            as.is = TRUE, row.names = 1)
        phenodata <- cbind(phenodata, pdata[rownames(phenodata), 
            ])
    }
    counts <- matrix(nrow = nrow(bins), ncol = length(bamnames), 
        dimnames = list(rownames(bins), bamnames))
    for (i in 1:length(bamfiles)) {
        counts[, i] <- .binReadCountsPerSample(bins, bamfiles[i], 
            path, cache, samtools, f, F, q, maxChunk)
        gc(FALSE)
    }
    phenodata$reads <- apply(counts, 2, sum)
    dat <- list(phenodata = phenodata, bins = bins, counts = counts)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

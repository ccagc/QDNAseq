\name{getBins}
\alias{getBins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getBins
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getBins(binsize, genome = "hg19", cache = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{binsize}{
%%     ~~Describe \code{binsize} here~~
}
  \item{genome}{
%%     ~~Describe \code{genome} here~~
}
  \item{cache}{
%%     ~~Describe \code{cache} here~~
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
function (binsize, genome = "hg19", cache = TRUE) 
{
    genome.build <- as.integer(gsub("[^0-9]", "", genome))
    if (genome.build \%in\% c(19, 37)) {
        genome.name <- "hg19"
    }
    else if (genome.build \%in\% c(18, 36)) {
        genome.name <- "hg18"
    }
    else {
        stop("Unknown genome: ", genome)
    }
    filename <- paste("QDNAseq.", genome.name, ".", binsize, 
        "kbp.txt.gz", sep = "")
    if (cache) {
        localfile <- file.path(path.expand("~"), ".QDNAseq", 
            filename)
    }
    else {
        localfile <- tempfile()
    }
    if (cache & !file.exists(file.path(path.expand("~"), ".QDNAseq"))) 
        dir.create(file.path(path.expand("~"), ".QDNAseq"))
    if (!file.exists(localfile)) {
        remotefile <- paste("http://cdn.bitbucket.org/ilarischeinin/qdnaseq/downloads/", 
            filename, sep = "")
        if (download.file(remotefile, localfile, quiet = TRUE) != 
            0) 
            stop("Annotations not found on server for genome ", 
                genome, " and bin size ", binsize, ". Please generate them first.")
    }
    bins <- read.table(localfile, header = TRUE, sep = "\t", 
        as.is = TRUE)
    rownames(bins) <- paste(bins$chromosome, ":", bins$start, 
        "-", bins$end, sep = "")
    bins
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

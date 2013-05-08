\name{poolReadCounts}
\alias{poolReadCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
poolReadCounts(obj, samples)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
%%     ~~Describe \code{obj} here~~
}
  \item{samples}{
%%     ~~Describe \code{samples} here~~
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
function (obj, samples) 
{
    phenodata <- obj[["phenodata"]]
    bins <- obj[["bins"]]
    counts <- obj[["counts"]]
    if ("corrected" \%in\% names(obj)) 
        corrected <- obj[["corrected"]]
    if (length(samples) != nrow(phenodata)) 
        stop("Vector samples must be of equal length as there are samples.")
    newsamples <- sort(unique(samples))
    oldsamples <- phenodata$name
    if (length(newsamples) == nrow(phenodata)) 
        return(obj)
    newphenodata <- data.frame(name = newsamples, reads = NA, 
        loess.span = NA, loess.family = NA, row.names = newsamples, 
        stringsAsFactors = FALSE)
    newcounts <- matrix(nrow = nrow(counts), ncol = length(newsamples), 
        dimnames = list(rownames(counts), newsamples))
    if ("corrected" \%in\% names(obj)) 
        newcorrected <- matrix(nrow = nrow(counts), ncol = length(newsamples), 
            dimnames = list(rownames(counts), newsamples))
    for (newsample in newsamples) {
        replicates <- samples == newsample
        newcounts[, x] <- rowSums(counts[, replicates, drop = FALSE])
        if ("corrected" \%in\% names(obj)) 
            newcorrected[, x] <- rowSums(corrected[, replicates, 
                drop = FALSE])
        oldpheno <- dat[["phenodata"]][replicates, ]
        newpheno[x, "reads"] <- sum(oldpheno$reads)
        if (length(unique(oldpheno$loess.span)) == 1) 
            newpheno[x, "loess.span"] <- oldpheno[1, "loess.span"]
        if (length(unique(oldpheno$loess.family)) == 1) 
            newpheno[x, "loess.family"] <- oldpheno[1, "loess.family"]
    }
    if ("corrected" \%in\% names(obj)) 
        return(list(phenodata = newphenodata, bins = bins, counts = newcounts, 
            corrected = newcorrected))
    list(phenodata = newphenodata, bins = bins, counts = newcounts)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

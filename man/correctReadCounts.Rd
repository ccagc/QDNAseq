\name{correctReadCounts}
\alias{correctReadCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{correctReadCounts
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
correctReadCounts(obj, span = 0.65, family = "symmetric", plotting = FALSE, condition = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
%%     ~~Describe \code{obj} here~~
}
  \item{span}{
%%     ~~Describe \code{span} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{plotting}{
%%     ~~Describe \code{plotting} here~~
}
  \item{condition}{
%%     ~~Describe \code{condition} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (obj, span = 0.65, family = "symmetric", plotting = FALSE, 
    condition = TRUE, ...) 
{
    phenodata <- obj[["phenodata"]]
    bins <- obj[["bins"]]
    counts <- obj[["counts"]]
    cat("Performing correction for GC content and mappability:\n")
    if (length(span) == 1) 
        span <- rep(span, ncol(counts))
    if (length(family) == 1) 
        family <- rep(family, ncol(counts))
    if (length(span) != ncol(counts)) 
        stop("span has to be either a single value or a vector the same length as the number of columns in counts.")
    if (length(family) != ncol(counts)) 
        stop("family has to be either a single value or a vector the same length as the number of columns in counts.")
    if (length(condition) == 1) {
        if (condition) {
            condition <- bins$chromosome \%in\% as.character(1:22) & 
                bins$chromosome == c(bins$chromosome[-1], "0")
            if ("blacklist" \%in\% colnames(bins)) 
                condition <- condition & bins$blacklist == 0
        }
        else {
            condition <- rep(TRUE, nrow(bins))
        }
    }
    if (length(condition) != nrow(bins)) 
        stop("condition has to be either a single value or vector the same length as the number of rows in bins.")
    used.span <- rep(NA, ncol(counts))
    used.family <- rep(NA, ncol(counts))
    corrected <- matrix(nrow = nrow(counts), ncol = ncol(counts), 
        dimnames = dimnames(counts))
    residuals <- matrix(nrow = nrow(counts), ncol = ncol(counts), 
        dimnames = dimnames(counts))
    gc <- round(bins$gc)
    mappability <- round(bins$mappability)
    median.counts <- aggregate(counts[condition], by = list(gc = gc[condition], 
        mappability = mappability[condition]), median)
    median.counts <- median.counts[!is.na(median.counts$gc), 
        ]
    rownames(median.counts) <- paste(median.counts$gc, "-", median.counts$mappability, 
        sep = "")
    for (i in 1:ncol(counts)) {
        if (is.na(span[i]) && is.na(family[i])) {
            cat("\tSkipping correction for sample ", colnames(counts)[i], 
                "...\n", sep = "")
            next
        }
        cat("\tUsing span=", span[i], "\tand family=", family[i], 
            ",\tcorrecting sample ", colnames(counts)[i], "...\n", 
            sep = "")
        vals <- median.counts[, i + 2]
        corvals <- counts[, i]
        try({
            l <- loess(vals ~ median.counts$gc * median.counts$mappability, 
                span = span[i], family = family[i])
            fit <- l$fitted
            names(fit) <- rownames(median.counts)
            residuals[, i] <- corvals - fit[paste(gc, "-", mappability, 
                sep = "")]
            correction <- median(fit, na.rm = TRUE) - fit
            corvals <- corvals + correction[paste(gc, "-", mappability, 
                sep = "")]
            corvals <- corvals - min(corvals, na.rm = TRUE)
            used.span[i] <- span[i]
            used.family[i] <- family[i]
        }, silent = TRUE)
        corrected[, i] <- corvals
    }
    corrected <- round(corrected, digits = 2)
    phenodata$loess.span <- used.span
    phenodata$loess.family <- used.family
    cat("Done.\n")
    list(phenodata = phenodata, bins = bins, counts = counts, 
        corrected = corrected, residuals = residuals)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

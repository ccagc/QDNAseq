%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  highlightFilters.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{highlightFilters}
\alias{highlightFilters}


\alias{highlightFilters,QDNAseqSignals-method}

\title{Highlights data points in a plotted profile to evaluate filtering}

\usage{
highlightFilters(object, col="red", residual=NA, blacklist=NA, mappability=NA, bases=NA,
  type="union", logTransform=TRUE, ...)
}

\description{
    Highlights data points in a plotted profile to evaluate filtering.
}

\arguments{
    \item{object}{A \code{\link{QDNAseqCopyNumbers}} object.}
    \item{col}{The color used for highlighting.}
    \item{residual}{Either a \code{\link[base]{logical}} specifying whether to filter based on
        loess residuals of the calibration set, or if a \code{\link[base]{numeric}}, the cutoff
        as number of standard deviations estimated with
        \code{\link[matrixStats]{madDiff}} to use for. Default is \code{\link[base:logical]{TRUE}}, which
        corresponds to 4.0 standard deviations.}
    \item{blacklist}{Either a \code{\link[base]{logical}} specifying whether to filter based on
        overlap with blacklisted regions, or if numeric, the maximum
        percentage of overlap allowed. Default is \code{\link[base:logical]{TRUE}}, which corresponds to
        no overlap allowd (i.e. value of 0).}
    \item{mappability}{A \code{\link[base]{numeric}} in \eqn{[0,100]} to specify filtering out
        bins with mappabilities lower than the number specified. NA (default)
        or \code{\link[base:logical]{FALSE}} will not filter based on mappability.}
    \item{bases}{A \code{\link[base]{numeric}} specifying the minimum percentage of characterized
        bases (not Ns) in the reference genome sequence. NA (default) or
        \code{\link[base:logical]{FALSE}} will not filted based on uncharacterized bases.}
    \item{type}{When specifying multiple filters (\code{residual},
        \code{blacklist}, \code{mappability}, \code{bases}), whether to
        highlight their \code{union} (default) or \code{intersection}.}
    \item{logTransform}{Whether the plotted data was log-transformed or not,
        should be the same value that was specified for \code{\link{plot}}().}
    \item{...}{Further arguments to \code{\link[graphics]{points}}.}
}

\examples{
data(LGG150)
readCounts <- LGG150
plot(readCounts)
highlightFilters(readCounts, residual=TRUE, blacklist=TRUE)
}

\author{Ilari Scheinin}


\keyword{aplot}

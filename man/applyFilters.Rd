%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  applyFilters.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{applyFilters}
\alias{applyFilters}


\alias{applyFilters,QDNAseqReadCounts-method}

\title{Adjusts the filtering on which bins are used}

\usage{
applyFilters(object, residual=TRUE, blacklist=TRUE, mappability=NA, bases=NA,
  filterAllosomes=TRUE)
}

\description{
    Adjusts the filtering on which bins are used.
}

\arguments{
    \item{object}{A \code{\link{QDNAseqReadCounts}} object.}
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
    \item{filterAllosomes}{A \code{\link[base]{logical}} specifying whether to filter out the sex
        chromosomes. Default is \code{\link[base:logical]{TRUE}}.}
}

\value{
    Returns a \code{\link{QDNAseqReadCounts}} object with updated filtering.
}

\examples{
data(LGG150)
readCounts <- LGG150
readCountsFiltered <- applyFilters(readCounts)
}

\author{Ilari Scheinin}


\keyword{manip}

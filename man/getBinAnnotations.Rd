%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  getBinAnnotations.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getBinAnnotations}
\alias{getBinAnnotations}


\title{Gets bin annotation data for a particular bin size}

\usage{
getBinAnnotations(binsize, genome="hg19", cache=TRUE, force=FALSE)
}

\description{
 Gets bin annotation data for a particular bin size.
}

\arguments{
  \item{binsize}{A \code{\link[base]{numeric}} scalar specifying the width of the bins
   in units of kbp (1000 base pairs), e.g. \code{binsize=15} corresponds
   to 15 kbp bins.}
  \item{genome}{A \code{\link[base]{character}} string specify the genome and genome version
   to be used.}
  \item{cache}{If \code{\link[base:logical]{TRUE}}, the retrieved bin annotation data is cached
   on the file system, otherwise not.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, the bin anonnation data is retrieved/calculated
   regardless of it already exists in the cache or not.}
}

\value{
  Returns a \code{\link[Biobase]{AnnotatedDataFrame}} object.
}

\author{Ilari Scheinin}

\seealso{
  \code{\link{createBins}}().
}


\keyword{IO}

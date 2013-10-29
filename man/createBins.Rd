%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  getBinAnnotations.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{createBins}
\alias{createBins}


\alias{calculateMappability}
\alias{calculateBlacklist}
\alias{iterateResiduals}

\title{Builds bin annotation data for a particular bin size}

\usage{
createBins(bsgenome, binSize, ignoreUnderscored=TRUE, ignoreMitochondria=TRUE)
}

\description{
 Builds bin annotation data for a particular bin size.
}

\arguments{
  \item{bsgenome}{A BSgenome ...}
  \item{binSize}{A \code{\link[base]{numeric}} scalar specifying the width of the bins
   in units of kbp (1000 base pairs), e.g. \code{binSize=15} corresponds
   to 15 kbp bins.}
  \item{ignoreUnderscored}{Whether to ignore sequences with underscores
    in their names ...}
  \item{ignoreMitochondria}{Wheter to ignore mitochondrial DNA  ...}
}

\value{
  Returns ...
}

\author{Ilari Scheinin}

\seealso{
  \code{\link{getBinAnnotations}}().
}


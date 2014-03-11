%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  plot-methods.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{noisePlot}
\alias{noisePlot}


\alias{noisePlot,QDNAseqReadCounts,missing-method}

\title{Plot noise as a function of sequence depth}

\usage{
noisePlot(x, y, ...)
}

\description{
 Plot noise as a function of sequence depth.
}

\arguments{
  \item{x}{A \code{\link{QDNAseqReadCounts}} object.}
  \item{y}{missing}
  \item{...}{Further arguments to \code{\link[graphics]{plot}} and
    \code{\link[graphics]{text}}.}
}

\examples{
  data(LGG150)
  readCounts <- LGG150
  readCountsFiltered <- applyFilters(readCounts)
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  noisePlot(readCountsFiltered)
}

\author{Ilari Scheinin}


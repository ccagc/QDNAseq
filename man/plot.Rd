%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  plot-methods.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plot}
\alias{plot}


\alias{plot,QDNAseqSignals,missing-method}

\title{Plot copy number profile}

\usage{
plot(x, y, ...)
}

\description{
 Plot copy number profile.
}

\arguments{
  \item{x}{A \code{\link{QDNAseqReadCounts}} or \code{\link{QDNAseqCopyNumbers}} object.}
  \item{y}{missing}
  \item{...}{...}
}

\examples{
  data(LGG150)
  readCounts <- LGG150
  readCountsFiltered <- applyFilters(readCounts)
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  copyNumbers <- correctBins(readCountsFiltered)
  plot(copyNumbers)
}

\author{Ilari Scheinin}


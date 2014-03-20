%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  plot-methods.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{frequencyPlot}
\alias{frequencyPlot}


\alias{frequencyPlot,QDNAseqCopyNumbers,missing-method}

\title{Plot copy number aberration frequencies}

\usage{
frequencyPlot(x, y, ...)
}

\description{
    Plot copy number aberration frequencies.
}

\arguments{
    \item{x}{A \code{\link{QDNAseqCopyNumbers}} object with \code{calls} data.}
    \item{y}{missing}
    \item{...}{...}
}

\examples{
data(LGG150)
readCounts <- LGG150
readCountsFiltered <- applyFilters(readCounts)
readCountsFiltered <- estimateCorrection(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
copyNumbersSegmented <- segmentBins(copyNumbersSmooth)
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
copyNumbersCalled <- callBins(copyNumbersSegmented)
frequencyPlot(copyNumbersCalled)
}

\author{Ilari Scheinin}


\keyword{hplot}

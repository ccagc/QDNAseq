%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  callBins.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{callBins}
\alias{callBins}


\alias{callBins,QDNAseqCopyNumbers-method}

\title{Call aberrations from segmented copy number data}

\usage{
callBins(object, ...)
}

\description{
 Call aberrations from segmented copy number data.
}

\arguments{
  \item{object}{An object of class QDNAseqCopyNumbers}
  \item{...}{Additional arguments passed to \code{\link[CGHcall]{CGHcall}}.}
}

\value{
  Returns an object of class QDNAseqCopyNumbers with segmentation results
    added.
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
}

\author{Ilari Scheinin}

\seealso{
  Internally, \code{\link[CGHcall]{CGHcall}} and \code{\link[CGHcall]{ExpandCGHcall}} of
    the \pkg{DNAcopy} package are used.
}


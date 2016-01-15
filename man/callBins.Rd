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
callBins(object, organism=c("human", "other"), method=c("CGHcall", "cutoff"),
  cutoffLoss=-0.8, cutoffGain=0.5, ...)
}

\description{
    Call aberrations from segmented copy number data.
}

\arguments{
    \item{object}{An object of class QDNAseqCopyNumbers}
    \item{organism}{Either \dQuote{human} or \dQuote{other}, see manual page
        for \code{\link[CGHcall]{CGHcall}} for more details. This is only used for
        chromosome arm information when \dQuote{prior} is set to \dQuote{all}
        or \dQuote{auto} (and samplesize > 20). Ignored when \code{method} is
        not \dQuote{CGHcall}.}
    \item{method}{Calling method to use. Options currently implemented are:
        \dQuote{CGHcall} or \dQuote{cutoff}.}
    \item{cutoffLoss}{When method=\dQuote{cutoff}, the (log2-transformed)
        threshold to use for losses. Default is -0.8. (The expected log2-ratio
        for a single copy loss present in all cells is -1.0.)}
    \item{cutoffGain}{When method=\dQuote{cutoff}, the (log2-transformed)
        threshold to use for gains. Default is 0.5. (The expected log2-ratio
        for a single copy gain present in all cells is 0.6.)}
    \item{...}{Additional arguments passed to \code{\link[CGHcall]{CGHcall}}.}
}

\details{
    Chromosomal aberrations are called with \pkg{CGHcall}. It has been
    developed for the analysis of series of cancer samples, and uses a model
    that contains both gains and losses. If used on a single sample, or
    especially only on a subset of chromosomes, or especially on a single
    non-cancer sample, it may fail.
}

\value{
    Returns an object of class \code{\link{QDNAseqCopyNumbers}} with calling
    results added.
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
    the \pkg{CGHcall} package are used.
}


\keyword{manip}

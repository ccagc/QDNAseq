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
  cutoffs=log2(c(deletion = 0.5, loss = 1.5, gain = 2.5, amplification = 10)/2), ...)
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
    \item{cutoffs}{When method=\dQuote{cutoff}, a numeric vector of
        (log2-transformed) thresholds to use for calling. At least one
        positive and one negative value must be provided. The smallest
        positive value is used as the cutoff for calling gains, and the
        negative value closest to zero is used as the cutoff for losses. If a
        second positive value is provided, it is used as the cutoff for
        amplifications. And if a second negative value is provided, it is used
        as the cutoff for homozygous deletions.}
    \item{...}{Additional arguments passed to \code{\link[CGHcall]{CGHcall}}.}
}

\details{
    By default, chromosomal aberrations are called with \pkg{CGHcall}. It has
    been developed for the analysis of series of cancer samples, and uses a
    model that contains both gains and losses. If used on a single sample, or
    especially only on a subset of chromosomes, or especially on a single
    non-cancer sample, it may fail, but method \dQuote{cutoff} can be used
    instead.

    When using method \dQuote{cutoff}, the default values assume a uniform
    cell population and correspond to thresholds of (assuming a diploid
    genome) 0.5, 1.5, 2.5, and 10 copies to distinguish between homozygous
    deletions, (hemizygous) losses, normal copy number, gains, and
    amplifications, respectively. When using with cancer samples, these values
    might require adjustments to account for tumor cell percentage.
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
    the \pkg{CGHcall} package are used when method=\dQuote{CGHcall}.
}


\keyword{manip}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  estimateCorrection.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{estimateCorrection}
\alias{estimateCorrection}


\alias{estimateCorrection,QDNAseqReadCounts-method}

\title{Estimate correction to read counts for GC content and mappability}

\usage{
estimateCorrection(object, span=0.65, family="symmetric", adjustIncompletes=TRUE,
  maxIter=1, cutoff=4, ...)
}

\description{
 Estimate correction to read counts for GC content and mappability.
}

\arguments{
  \item{object}{An \code{\link{QDNAseqReadCounts}} object with \code{counts} data.}
  \item{span}{For \code{\link[stats]{loess}}, the parameter alpha which controls the
    degree of smoothing.}
  \item{family}{For \code{\link[stats]{loess}}, if "gaussian" fitting is by
    least-squares, and if "symmetric" a re-descending M estimator is
    used with Tukey's biweight function.}
  \item{adjustIncompletes}{A boolean(1) specifying whether \code{counts} for
    bins with uncharacterized nucleotides (N's) in their reference genome
    sequence should be adjusted by dividing them with the percentage of
    characterized (A, C, G, T) nucleotides. Defaults to \code{\link[base:logical]{TRUE}}.}
  \item{maxIter}{An integer(1) specifying the maximum number of iterations to perform, default is 1. If larger, after the first loess fit, bins with median residuals larger than \code{cutoff} are removed, and the fitting repeated until the list of bins to use stabilizes or after \code{maxIter} iterations.}
  \item{cutoff}{A numeric(1) specifying the number of standard deviations (as
    estimated with \code{\link[matrixStats]{madDiff}}) the cutoff for removal of bins
    with median residuals larger than the cutoff. Not used if \code{maxIter=1}
    (default).}
  \item{...}{Additional aguments passed to \code{\link[stats]{loess}}.}
}

\value{
  Returns a \code{\link{QDNAseqReadCounts}} object with the assay data element
  \code{fit} added.
}

\examples{
  data(LGG150)
  readCounts <- LGG150
  readCountsFiltered <- applyFilters(readCounts)
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
}

\author{Ilari Scheinin}

\seealso{
  Internally, \code{\link[stats]{loess}} is used to fit the regression model.
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  correctBins.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{correctBins}
\alias{correctBins}


\alias{correctBins,QDNAseqReadCounts-method}

\title{Correct binned read counts for GC content and mappability}

\usage{
correctBins(object, fit=NULL, method="ratio", adjustIncompletes=TRUE, ...)
}

\description{
 Correct binned read counts for GC content and mappability.
}

\arguments{
  \item{object}{An \code{\link{QDNAseqReadCounts}} object with \code{counts} data.}
  \item{fit}{...}
  \item{method}{...}
  \item{adjustIncompletes}{...}
  \item{...}{Additional arguments passed to \code{\link{estimateCorrection}}().}
}

\value{
  Returns a \code{\link{QDNAseqReadCounts}} object with the assay data element
  \code{corrected} added.
}

\author{Ilari Scheinin}

\seealso{
  Internally, \code{\link[stats]{loess}} is used to fit the regression model.
}


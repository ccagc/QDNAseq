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
  \item{span}{...}
  \item{family}{...}
  \item{adjustIncompletes}{...}
  \item{maxIter}{...}
  \item{cutoff}{...}
  \item{...}{Additional aguments passed to \code{\link[stats]{loess}}.}
}

\value{
  Returns a \code{\link{QDNAseqReadCounts}} object with the assay data element
  \code{fit} added.
}

\author{Ilari Scheinin}

\seealso{
  Internally, \code{\link[stats]{loess}} is used to fit the regression model.
}


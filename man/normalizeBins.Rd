%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeBins.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{normalizeBins}
\alias{normalizeBins}


\alias{normalizeBins,QDNAseqReadCounts-method}

\title{Normalizes binned read counts}

\usage{
normalizeBins(object, method="median", smoothOutliers=TRUE, logTransform=TRUE,
  logOffset=2^-10, force=FALSE, ...)
}

\description{
 Normalizes binned read counts.
}

\arguments{
  \item{object}{...}
  \item{method}{A \code{\link[base]{character}} string specifying ...}
  \item{smoothOutliers}{If \code{\link[base:logical]{TRUE}}, ..., otherwise, ...}
  \item{logTransform}{If \code{\link[base:logical]{TRUE}}, ..., otherwise, ...}
  \item{logOffset}{...}
  \item{force}{...}
  \item{...}{Not used.}
}

\value{
  Returns a named \code{\link[base]{list}} containing elements ...
}

\author{Ilari Scheinin}

\seealso{
  Internally, ...
}


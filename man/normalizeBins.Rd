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
  \item{object}{A QDNAseqReadCounts object.}
  \item{method}{A \code{\link[base]{character}} string specifying the normalization method.
    Choices are "none", "mean", "median" (default), or "mode". A partial
    string sufficient to uniquely identify the choice is permitted.}
  \item{smoothOutliers}{If \code{\link[base:logical]{TRUE}} (default), smooth.CNA from package DNAcopy is
    used detect and smooth outliers.}
  \item{logTransform}{If \code{\link[base:logical]{TRUE}} (default), data will be log2-transformed.}
  \item{logOffset}{An offset to be added prior to log2-transformation to
    avoid non-positive numbers. Ignored if logTransform=\code{\link[base:logical]{FALSE}}.}
  \item{force}{Running this function will remove possible segmentation and
    calling results. When they are present, running requires specifying
    force=\code{\link[base:logical]{TRUE}}.}
  \item{...}{Further arguments to DNAcopy::smooth.CNA.}
}

\value{
  Returns a QDNAseqReadCounts object with normalized data.
}

\author{Ilari Scheinin}


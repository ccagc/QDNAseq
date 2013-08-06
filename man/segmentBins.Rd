%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  segmentBins.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{segmentBins}
\alias{segmentBins}


\alias{segmentBins,QDNAseqReadCounts-method}

\title{Segments and calls total copy numbers}

\usage{
segmentBins(object, weights=FALSE, normalize=TRUE, inter=c(-0.1, 0.1), force=FALSE, ...)
}

\description{
 Segments and calls total copy numbers.
}

\arguments{
  \item{object}{...}
  \item{weights}{Either \code{\link[base:logical]{TRUE}} or a vector of weights. If \code{\link[base:logical]{TRUE}},
    loess residuals are used as weights.}
  \item{normalize}{...}
  \item{inter}{...}
  \item{force}{...}
  \item{...}{Additional arguments passed to \code{\link[DNAcopy]{segment}}.}
}

\value{
  Returns ...
}

\author{Ilari Scheinin}

\seealso{
  Internally, \code{\link[DNAcopy]{segment}} of the \pkg{DNAcopy} package,
  which implements the CBS method, is used to segment the data.
}


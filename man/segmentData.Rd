%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  segmentData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{segmentData}
\alias{segmentData}


\alias{segmentData,QDNAseqReadCounts,logical-method}
\alias{segmentData,cghRaw,missing-method}

\title{Segments and calls total copy numbers}

\usage{
segmentData(object, weights, ...)
}

\description{
 Segments and calls total copy numbers.
}

\arguments{
  \item{object}{...}
  \item{weights}{Either \code{\link[base:logical]{TRUE}} or a vector of weights. If \code{\link[base:logical]{TRUE}},
    loess residuals are used as weights.}
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

